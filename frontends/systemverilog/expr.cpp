/*! \file expr.cpp
 *  \brief SystemVerilogEncoder's expression-to-term conversion: the main
 *         expr_to_term() switch covering literals, operators, selects,
 *         conversions, and system-function calls.
 */
#include <string>

#include "frontends/systemverilog/ast_helpers.h"
#include "frontends/systemverilog/encoder.h"
#include "slang/ast/EvalContext.h"
#include "slang/ast/Expression.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/expressions/AssignmentExpressions.h"
#include "slang/ast/expressions/CallExpression.h"
#include "slang/ast/expressions/ConversionExpression.h"
#include "slang/ast/expressions/LiteralExpressions.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/expressions/Operator.h"
#include "slang/ast/expressions/OperatorExpressions.h"
#include "slang/ast/expressions/SelectExpressions.h"
#include "slang/ast/symbols/MemberSymbols.h"
#include "slang/ast/symbols/VariableSymbols.h"
#include "slang/ast/types/AllTypes.h"
#include "slang/ast/types/Type.h"
#include "slang/numeric/SVInt.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

Term SystemVerilogEncoder::expr_to_term(const slang::ast::Expression & expr)
{
  using namespace slang::ast;

  switch (expr.kind) {
    case ExpressionKind::NamedValue: {
      auto & nv = expr.as<NamedValueExpression>();
      return lookup_symbol(&canonicalize_modport_port(nv.symbol));
    }

    case ExpressionKind::HierarchicalValue: {
      // Cross-scope dotted read (e.g. `child_inst.reg`).  Slang has
      // already resolved the dotted path to the target ValueSymbol;
      // lookup_symbol finds its term in the appropriate scope
      // provided the referenced instance has been encoded already.
      // A modport-qualified access resolves to a ModportPortSymbol
      // proxy rather than the real symbol directly -- see
      // canonicalize_modport_port().
      auto & hv = expr.as<HierarchicalValueExpression>();
      return lookup_symbol(&canonicalize_modport_port(hv.symbol));
    }

    case ExpressionKind::LValueReference: {
      // Implicit self-reference produced by compound assignments
      // (`x &= y` -> `x = LValueReference & y`).  The owning
      // assignment handler must have stashed the current LHS term
      // before recursing into the RHS.
      if (!current_lvalue_term_) {
        throw PonoException(
            "SystemVerilogEncoder: LValueReference outside compound "
            "assignment");
      }
      return current_lvalue_term_;
    }

    case ExpressionKind::IntegerLiteral: {
      auto & lit = expr.as<IntegerLiteral>();
      uint64_t width = expr.type->getBitWidth();
      if (width == 0) width = 32;  // Default integer width.
      Sort sort = solver_->make_sort(BV, width);
      // Reinterpret the value as unsigned so that toString emits the raw
      // two's-complement bit pattern as a positive decimal.  Without
      // setSigned(false), signed-negative values would stringify as
      // "-N", which smt-switch's base-10 parser rejects.
      auto val = lit.getValue();
      val.setSigned(false);
      string val_str =
          val.toString(slang::LiteralBase::Decimal, /*includeBase=*/false);
      return solver_->make_term(val_str, sort, 10);
    }

    case ExpressionKind::UnbasedUnsizedIntegerLiteral: {
      auto & lit = expr.as<UnbasedUnsizedIntegerLiteral>();
      uint64_t width = expr.type->getBitWidth();
      if (width == 0) width = 1;
      Sort sort = solver_->make_sort(BV, width);
      auto val = lit.getValue();
      val.setSigned(false);
      string val_str =
          val.toString(slang::LiteralBase::Decimal, /*includeBase=*/false);
      return solver_->make_term(val_str, sort, 10);
    }

    case ExpressionKind::BinaryOp: {
      auto & binop = expr.as<BinaryExpression>();

      if (binop.op == BinaryOperator::WildcardEquality
          || binop.op == BinaryOperator::WildcardInequality) {
        // Must special-case *before* the generic eager left/right
        // conversion below: a right operand with unknown bits (e.g.
        // `4'b10??`) would otherwise hit the generic (wildcard-
        // unaware) IntegerLiteral case first and crash trying to hand
        // an X-containing decimal string to the solver, the same bug
        // casex/casez's item patterns had. Per the LRM, only the
        // *right* operand's X/Z bits are wildcards, so left is always
        // safe to convert normally; only convert right if we end up
        // needing it for the plain-equality fallback below. Reuses
        // the same (mask, value) technique as casex/casez in
        // process_statement()'s Case handling: build a mask with a 0
        // at each wildcard bit and compare (left & mask) == value,
        // ignoring exactly those positions. Falls back to plain
        // equality if the right operand isn't a usable compile-time
        // constant of 1-64 bits (mask/value are tracked in a
        // uint64_t) -- nothing else to wildcard against, and this
        // encoder's BV model has no way for a non-literal term to
        // hold an unknown bit at all.
        Term left = expr_to_term(binop.left());
        uint64_t result_width = expr.type->getBitWidth();
        Term eq;
        auto rhs_cv = binop.right().eval(eval_ctx());
        uint64_t rhs_w = binop.right().type->getBitWidth();
        if (!rhs_cv.bad() && rhs_cv.isInteger() && rhs_w > 0 && rhs_w <= 64) {
          auto & sv = rhs_cv.integer();
          uint64_t mask = 0, value = 0;
          for (uint64_t i = 0; i < rhs_w; ++i) {
            slang::logic_t bit = sv[static_cast<int32_t>(i)];
            if (!bit.isUnknown()) {
              mask |= (uint64_t{ 1 } << i);
              if (bit.value == 1) value |= (uint64_t{ 1 } << i);
            }
          }
          Sort rhs_sort = solver_->make_sort(BV, rhs_w);
          // Raw bit-pattern masks, not numeric values -- always
          // zero-extend.
          Term mask_term =
              resize_to(solver_->make_term(std::to_string(mask), rhs_sort, 10),
                        left->get_sort()->get_width(),
                        false);
          Term value_term =
              resize_to(solver_->make_term(std::to_string(value), rhs_sort, 10),
                        left->get_sort()->get_width(),
                        false);
          eq = solver_->make_term(
              Equal, solver_->make_term(BVAnd, left, mask_term), value_term);
        } else {
          Term right = expr_to_term(binop.right());
          right = resize_to(right,
                            left->get_sort()->get_width(),
                            binop.right().type->isSigned());
          eq = solver_->make_term(Equal, left, right);
        }
        Sort bv1 = solver_->make_sort(BV, 1);
        bool want_eq = binop.op == BinaryOperator::WildcardEquality;
        Term result =
            solver_->make_term(Ite,
                               eq,
                               solver_->make_term(want_eq ? 1 : 0, bv1),
                               solver_->make_term(want_eq ? 0 : 1, bv1));
        // Result is a plain 0/1 boolean -- always zero-extend.
        return resize_to(result, result_width, false);
      }

      Term left = expr_to_term(binop.left());
      Term right = expr_to_term(binop.right());

      // Both operands are signed iff the *whole* operation is signed
      // (SystemVerilog's usual arithmetic conversion rule: mixing a
      // signed operand with an unsigned one makes the whole operation
      // unsigned) -- checking both, rather than trusting they already
      // agree, doesn't rely on slang having inserted a unifying
      // conversion on every operator variant that reaches here.
      bool op_signed =
          binop.left().type->isSigned() && binop.right().type->isSigned();

      // Ensure operands have the same width, sign-extending rather
      // than zero-extending a narrower signed operand -- otherwise a
      // negative value becomes a large positive one before the
      // operator below ever sees it.
      uint64_t lw = left->get_sort()->get_width();
      uint64_t rw = right->get_sort()->get_width();
      uint64_t max_w = max(lw, rw);
      left = resize_to(left, max_w, op_signed);
      right = resize_to(right, max_w, op_signed);

      uint64_t result_width = expr.type->getBitWidth();
      Term result;

      switch (binop.op) {
        case BinaryOperator::Add:
          result = solver_->make_term(BVAdd, left, right);
          break;
        case BinaryOperator::Subtract:
          result = solver_->make_term(BVSub, left, right);
          break;
        case BinaryOperator::Multiply:
          result = solver_->make_term(BVMul, left, right);
          break;
        case BinaryOperator::Divide:
          result = solver_->make_term(op_signed ? BVSdiv : BVUdiv, left, right);
          break;
        case BinaryOperator::Mod:
          result = solver_->make_term(op_signed ? BVSrem : BVUrem, left, right);
          break;
        case BinaryOperator::BinaryAnd:
          result = solver_->make_term(BVAnd, left, right);
          break;
        case BinaryOperator::BinaryOr:
          result = solver_->make_term(BVOr, left, right);
          break;
        case BinaryOperator::BinaryXor:
          result = solver_->make_term(BVXor, left, right);
          break;
        case BinaryOperator::BinaryXnor: {
          Term xor_t = solver_->make_term(BVXor, left, right);
          result = solver_->make_term(BVNot, xor_t);
          break;
        }
        case BinaryOperator::Equality: {
          Term eq = solver_->make_term(Equal, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, eq, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::Inequality: {
          Term eq = solver_->make_term(Equal, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, eq, solver_->make_term(0, bv1), solver_->make_term(1, bv1));
          break;
        }
        case BinaryOperator::LessThan: {
          Term lt = solver_->make_term(op_signed ? BVSlt : BVUlt, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, lt, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::LessThanEqual: {
          Term le = solver_->make_term(op_signed ? BVSle : BVUle, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, le, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::GreaterThan: {
          Term gt = solver_->make_term(op_signed ? BVSlt : BVUlt, right, left);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, gt, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::GreaterThanEqual: {
          Term ge = solver_->make_term(op_signed ? BVSle : BVUle, right, left);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, ge, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::LogicalAnd: {
          // Logical AND: both operands nonzero.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term l_nz = solver_->make_term(
              Distinct, left, solver_->make_term(0, left->get_sort()));
          Term r_nz = solver_->make_term(
              Distinct, right, solver_->make_term(0, right->get_sort()));
          Term both = solver_->make_term(And, l_nz, r_nz);
          result = solver_->make_term(Ite,
                                      both,
                                      solver_->make_term(1, bv1),
                                      solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::LogicalOr: {
          Sort bv1 = solver_->make_sort(BV, 1);
          Term l_nz = solver_->make_term(
              Distinct, left, solver_->make_term(0, left->get_sort()));
          Term r_nz = solver_->make_term(
              Distinct, right, solver_->make_term(0, right->get_sort()));
          Term either = solver_->make_term(Or, l_nz, r_nz);
          result = solver_->make_term(Ite,
                                      either,
                                      solver_->make_term(1, bv1),
                                      solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::LogicalShiftLeft:
          result = solver_->make_term(BVShl, left, right);
          break;
        case BinaryOperator::LogicalShiftRight:
          result = solver_->make_term(BVLshr, left, right);
          break;
        case BinaryOperator::ArithmeticShiftLeft:
          result = solver_->make_term(BVShl, left, right);
          break;
        case BinaryOperator::ArithmeticShiftRight:
          result = solver_->make_term(BVAshr, left, right);
          break;
        case BinaryOperator::Power: {
          // Scoped to a compile-time-constant exponent (the
          // overwhelming majority of real synthesizable use, e.g.
          // `x**2`, `2**WIDTH` in a parameter expression) -- unrolled
          // into repeated multiplication at the operands' common
          // width, matching the truncating-wraparound semantics the
          // rest of this encoder already uses for arithmetic (the
          // uniform resize_to(result, result_width) below then
          // applies exactly as it does for every other operator). A
          // non-constant exponent would need real BV exponentiation,
          // which isn't part of the SMT BV theory and isn't worth a
          // barrel-multiplier-style encoding for how rarely it's used
          // with a runtime exponent.
          auto exp_cv = binop.right().eval(eval_ctx());
          if (exp_cv.bad() || !exp_cv.isInteger()) {
            throw PonoException(
                "SystemVerilogEncoder: '**' is only supported with a "
                "compile-time-constant exponent");
          }
          auto exp_opt = exp_cv.integer().as<uint32_t>();
          if (!exp_opt) {
            throw PonoException(
                "SystemVerilogEncoder: '**' exponent out of range");
          }
          result = solver_->make_term(1, left->get_sort());
          for (uint32_t i = 0; i < *exp_opt; ++i) {
            result = solver_->make_term(BVMul, result, left);
          }
          break;
        }
        case BinaryOperator::CaseEquality: {
          // No X/Z representation in this encoder's pure-BV model, so
          // case equality can never actually differ from logical
          // equality.
          Term eq = solver_->make_term(Equal, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, eq, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::CaseInequality: {
          Term eq = solver_->make_term(Equal, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, eq, solver_->make_term(0, bv1), solver_->make_term(1, bv1));
          break;
        }
        default:
          throw PonoException(
              "SystemVerilogEncoder: unsupported binary operator "
              + to_string(static_cast<int>(binop.op)));
      }
      // If the surrounding (context-determined) width is wider than
      // the operands' own common width, extend the already-computed
      // result rather than the operands -- self-determined operations
      // (e.g. wraparound on overflow) happen at the operands' natural
      // width first, exactly as the LRM specifies, and only the final
      // value is widened to fit the context.
      return resize_to(result, result_width, expr.type->isSigned());
    }

    case ExpressionKind::UnaryOp: {
      auto & unop = expr.as<UnaryExpression>();
      Term operand = expr_to_term(unop.operand());
      uint64_t result_width = expr.type->getBitWidth();

      Term result;
      switch (unop.op) {
        case UnaryOperator::BitwiseNot:
          result = solver_->make_term(BVNot, operand);
          break;
        case UnaryOperator::LogicalNot: {
          Sort bv1 = solver_->make_sort(BV, 1);
          Term is_zero = solver_->make_term(
              Equal, operand, solver_->make_term(0, operand->get_sort()));
          result = solver_->make_term(Ite,
                                      is_zero,
                                      solver_->make_term(1, bv1),
                                      solver_->make_term(0, bv1));
          break;
        }
        case UnaryOperator::Minus:
          result = solver_->make_term(BVNeg, operand);
          break;
        case UnaryOperator::BitwiseAnd: {
          // Reduction AND: result is 1 if all bits are 1.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term all_ones = solver_->make_term(
              Equal,
              operand,
              solver_->make_term(string(operand->get_sort()->get_width(), '1'),
                                 operand->get_sort(),
                                 2));
          result = solver_->make_term(Ite,
                                      all_ones,
                                      solver_->make_term(1, bv1),
                                      solver_->make_term(0, bv1));
          break;
        }
        case UnaryOperator::BitwiseOr: {
          // Reduction OR: result is 1 if any bit is 1.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term any_one = solver_->make_term(
              Distinct, operand, solver_->make_term(0, operand->get_sort()));
          result = solver_->make_term(Ite,
                                      any_one,
                                      solver_->make_term(1, bv1),
                                      solver_->make_term(0, bv1));
          break;
        }
        case UnaryOperator::BitwiseXor: {
          // Reduction XOR: parity of bits. For a BV of width n,
          // XOR all bits together.
          uint64_t w = operand->get_sort()->get_width();
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(Op(Extract, 0, 0), operand);
          for (uint64_t i = 1; i < w; i++) {
            Term bit = solver_->make_term(Op(Extract, i, i), operand);
            result = solver_->make_term(BVXor, result, bit);
          }
          break;
        }
        case UnaryOperator::Plus:
          // Unary `+` is a no-op per the LRM.
          result = operand;
          break;
        case UnaryOperator::BitwiseNand: {
          // Reduction NAND: NOT(AND-reduce) -- same all-ones check as
          // BitwiseAnd above, with the Ite branches swapped.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term all_ones = solver_->make_term(
              Equal,
              operand,
              solver_->make_term(string(operand->get_sort()->get_width(), '1'),
                                 operand->get_sort(),
                                 2));
          result = solver_->make_term(Ite,
                                      all_ones,
                                      solver_->make_term(0, bv1),
                                      solver_->make_term(1, bv1));
          break;
        }
        case UnaryOperator::BitwiseNor: {
          // Reduction NOR: NOT(OR-reduce) -- same any-one check as
          // BitwiseOr above, with the Ite branches swapped.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term any_one = solver_->make_term(
              Distinct, operand, solver_->make_term(0, operand->get_sort()));
          result = solver_->make_term(Ite,
                                      any_one,
                                      solver_->make_term(0, bv1),
                                      solver_->make_term(1, bv1));
          break;
        }
        case UnaryOperator::BitwiseXnor: {
          // Reduction XNOR: NOT(XOR-reduce parity) -- same bit-by-bit
          // XOR fold as BitwiseXor above, negated at the end.
          uint64_t w = operand->get_sort()->get_width();
          Term parity = solver_->make_term(Op(Extract, 0, 0), operand);
          for (uint64_t i = 1; i < w; i++) {
            Term bit = solver_->make_term(Op(Extract, i, i), operand);
            parity = solver_->make_term(BVXor, parity, bit);
          }
          result = solver_->make_term(BVNot, parity);
          break;
        }
        default:
          throw PonoException(
              "SystemVerilogEncoder: unsupported unary operator "
              + to_string(static_cast<int>(unop.op)));
      }
      // See the analogous BinaryOp comment above: widen the computed
      // result to fit a wider context-determined width, sign-extending
      // if the result is itself signed.
      return resize_to(result, result_width, expr.type->isSigned());
    }

    case ExpressionKind::Streaming: {
      // Scoped to a single stream with no `with` sub-range clause --
      // covers real synthesizable usage of streaming concatenation
      // (reversing/regrouping one packed value's bits or byte-lanes,
      // e.g. `{<<{a}}`), not the LRM's full generality (multiple
      // streams, `with` ranges, dynamically-sized queues/arrays).
      auto & sc = expr.as<StreamingConcatenationExpression>();
      if (sc.streams().size() != 1 || sc.streams()[0].withExpr) {
        throw PonoException(
            "SystemVerilogEncoder: unsupported streaming concatenation "
            "shape (only a single stream with no `with` range is "
            "supported)");
      }
      Term base = expr_to_term(*sc.streams()[0].operand);
      uint64_t slice = sc.getSliceSize();
      if (slice == 0) {
        // `{>>{x}}` (left-to-right): a single stream is unchanged.
        return base;
      }
      uint64_t w = base->get_sort()->get_width();
      if (w % slice != 0) {
        throw PonoException(
            "SystemVerilogEncoder: streaming concatenation slice size "
            "does not evenly divide operand width");
      }
      // `{<<{x}}` (right-to-left): reassemble slice-wide chunks in
      // reverse order, taken from x's LSB to MSB.
      Term result;
      for (uint64_t lo = 0; lo < w; lo += slice) {
        Term piece = solver_->make_term(Op(Extract, lo + slice - 1, lo), base);
        result = result ? solver_->make_term(Concat, result, piece) : piece;
      }
      return result;
    }

    case ExpressionKind::Conversion: {
      // Widening a *signed* source value must sign-extend (replicate
      // the top bit) rather than zero-extend, or a negative value
      // becomes a large positive one in the wider representation.
      // The source operand's own signedness decides this, not the
      // target type's: converting an unsigned value into a wider
      // signed type still zero-extends, since there is no sign bit in
      // the source to replicate.
      auto & conv = expr.as<ConversionExpression>();
      Term inner = expr_to_term(conv.operand());
      uint64_t target_width = expr.type->getBitWidth();
      return resize_to(inner, target_width, conv.operand().type->isSigned());
    }

    case ExpressionKind::Concatenation: {
      auto & concat = expr.as<ConcatenationExpression>();
      auto operands = concat.operands();
      if (operands.empty()) {
        throw PonoException("SystemVerilogEncoder: empty concatenation");
      }
      // Concatenate from MSB (first) to LSB (last).
      Term result = expr_to_term(*operands[0]);
      for (size_t i = 1; i < operands.size(); i++) {
        Term next = expr_to_term(*operands[i]);
        result = solver_->make_term(Concat, result, next);
      }
      return result;
    }

    case ExpressionKind::Replication: {
      auto & repl = expr.as<ReplicationExpression>();
      Term inner = expr_to_term(repl.concat());
      // The count should be a compile-time constant.
      auto count_cv = repl.count().getConstant();
      if (!count_cv) {
        throw PonoException(
            "SystemVerilogEncoder: non-constant replication count");
      }
      auto count_opt = count_cv->integer().as<uint32_t>();
      if (!count_opt) {
        throw PonoException("SystemVerilogEncoder: invalid replication count");
      }
      uint32_t count = *count_opt;
      if (count == 0) {
        throw PonoException("SystemVerilogEncoder: zero replication count");
      }
      Term result = inner;
      for (uint32_t i = 1; i < count; i++) {
        result = solver_->make_term(Concat, result, inner);
      }
      return result;
    }

    case ExpressionKind::ElementSelect: {
      // Element select: `a[i]`.  For packed arrays the element width
      // can be more than one bit; both the constant-index case (Extract
      // over the full slice) and the dynamic case (shift by idx*elem_w
      // then extract elem_w bits) handle arbitrary element widths.
      auto & sel = expr.as<ElementSelectExpression>();
      Term val = expr_to_term(sel.value());
      auto & sel_expr = sel.selector();
      uint64_t elem_w = expr.type->getBitWidth();
      if (elem_w == 0) elem_w = 1;

      // Try to evaluate the selector as a constant -- including the
      // case where the index is a loop counter bound via eval_ctx.
      std::optional<uint64_t> idx_const;
      if (sel_expr.getConstant()) {
        idx_const = sel_expr.getConstant()->integer().as<uint64_t>();
      } else {
        auto cv = sel_expr.eval(eval_ctx());
        if (cv.isInteger()) idx_const = cv.integer().as<uint64_t>();
      }

      if (idx_const) {
        uint64_t low = *idx_const * elem_w;
        uint64_t high = low + elem_w - 1;
        return solver_->make_term(Op(Extract, high, low), val);
      }

      // Dynamic select: shift right by (idx * elem_w) bits, then
      // extract the bottom `elem_w` bits.
      uint64_t val_w = val->get_sort()->get_width();
      Sort val_sort = solver_->make_sort(BV, val_w);
      // Index arithmetic for a shift amount, not a value -- zero-
      // extend regardless of the index expression's own signedness.
      Term idx = expr_to_term(sel_expr);
      idx = resize_to(idx, val_w, false);
      Term shift_amount = idx;
      if (elem_w != 1) {
        Term elem_w_term = solver_->make_term(elem_w, val_sort);
        shift_amount = solver_->make_term(BVMul, idx, elem_w_term);
      }
      Term shifted = solver_->make_term(BVLshr, val, shift_amount);
      return solver_->make_term(Op(Extract, elem_w - 1, 0), shifted);
    }

    case ExpressionKind::RangeSelect: {
      // Range select: a[hi:lo]
      auto & sel = expr.as<RangeSelectExpression>();
      Term val = expr_to_term(sel.value());
      auto & left_expr = sel.left();
      auto & right_expr = sel.right();

      // Both bounds should be compile-time constants for synthesizable code.
      if (!left_expr.getConstant() || !right_expr.getConstant()) {
        throw PonoException(
            "SystemVerilogEncoder: non-constant range select bounds");
      }
      auto hi_opt = left_expr.getConstant()->integer().as<uint64_t>();
      auto lo_opt = right_expr.getConstant()->integer().as<uint64_t>();
      if (!hi_opt || !lo_opt) {
        throw PonoException(
            "SystemVerilogEncoder: invalid range select bounds");
      }
      uint64_t hi = *hi_opt;
      uint64_t lo = *lo_opt;
      if (hi < lo) swap(hi, lo);
      return solver_->make_term(Op(Extract, hi, lo), val);
    }

    case ExpressionKind::MemberAccess: {
      // Packed-struct field read (`s.field`): extract the field's bit
      // range, using the FieldSymbol's own bitOffset (LSB-relative,
      // per packed-struct layout) and declared width.
      auto & ma = expr.as<MemberAccessExpression>();
      if (ma.member.kind != SymbolKind::Field) {
        throw PonoException(
            "SystemVerilogEncoder: unsupported member access on "
            + std::string(ma.member.name));
      }
      Term base = expr_to_term(ma.value());
      auto & field = ma.member.as<FieldSymbol>();
      uint64_t w = field.getType().getBitWidth();
      uint64_t lo = field.bitOffset;
      uint64_t hi = lo + w - 1;
      return solver_->make_term(Op(Extract, hi, lo), base);
    }

    case ExpressionKind::StructuredAssignmentPattern: {
      // Packed-struct construction (`'{field: value, ...}`): every
      // field must be given a named setter -- concatenate them MSB
      // first (declaration order), matching packed-struct layout.
      auto & pat = expr.as<StructuredAssignmentPatternExpression>();
      auto & canon = expr.type->getCanonicalType();
      if (canon.kind != SymbolKind::PackedStructType) {
        throw PonoException(
            "SystemVerilogEncoder: unsupported assignment pattern target "
            "type");
      }
      auto & st = canon.as<PackedStructType>();
      std::unordered_map<const Symbol *, const Expression *> setters;
      for (auto & ms : pat.memberSetters) {
        setters[&*ms.member] = &*ms.expr;
      }
      std::vector<Term> parts;
      for (auto & m : st.members()) {
        if (m.kind != SymbolKind::Field) continue;
        auto & field = m.as<FieldSymbol>();
        auto sit = setters.find(&m);
        if (sit == setters.end()) {
          throw PonoException(
              "SystemVerilogEncoder: assignment pattern missing field '"
              + std::string(field.name) + "'");
        }
        Term val = resize_to(expr_to_term(*sit->second),
                             field.getType().getBitWidth(),
                             sit->second->type->isSigned());
        parts.push_back(val);
      }
      if (parts.empty()) {
        throw PonoException(
            "SystemVerilogEncoder: assignment pattern for empty struct");
      }
      Term result = parts[0];
      for (size_t i = 1; i < parts.size(); ++i) {
        result = solver_->make_term(Concat, result, parts[i]);
      }
      return result;
    }

    case ExpressionKind::ConditionalOp: {
      // `cond1 &&& cond2 &&& ... ? a : b` (LRM 11.4.11): more than one
      // `&&&`-joined condition is legal outside case/if context too.
      // A `matches` pattern on any condition introduces destructuring
      // bind semantics this encoder doesn't implement -- throw rather
      // than silently evaluate only the plain boolean part of it.
      auto & ternary = expr.as<ConditionalExpression>();
      Term bool_cond;
      for (auto & c : ternary.conditions) {
        if (c.pattern) {
          throw PonoException(
              "SystemVerilogEncoder: pattern-matching ternary condition "
              "('... matches ...') is not supported");
        }
        Term c_term = expr_to_term(*c.expr);
        uint64_t cw = c_term->get_sort()->get_width();
        Term c_bool =
            (cw == 1)
                ? solver_->make_term(
                      Equal,
                      c_term,
                      solver_->make_term(1, solver_->make_sort(BV, 1)))
                : solver_->make_term(Distinct,
                                     c_term,
                                     solver_->make_term(0, c_term->get_sort()));
        bool_cond =
            bool_cond ? solver_->make_term(And, bool_cond, c_bool) : c_bool;
      }
      Term then_val = expr_to_term(ternary.left());
      Term else_val = expr_to_term(ternary.right());

      // Ensure then/else have the same width, sign-extending rather
      // than zero-extending a narrower signed branch (see BinaryOp's
      // op_signed handling above for why).
      bool branches_signed =
          ternary.left().type->isSigned() && ternary.right().type->isSigned();
      uint64_t tw = then_val->get_sort()->get_width();
      uint64_t ew = else_val->get_sort()->get_width();
      uint64_t max_w = max(tw, ew);
      then_val = resize_to(then_val, max_w, branches_signed);
      else_val = resize_to(else_val, max_w, branches_signed);

      return solver_->make_term(Ite, bool_cond, then_val, else_val);
    }

    case ExpressionKind::Call: {
      // `$signed`/`$unsigned` are a pure type reinterpretation;
      // `$past`/`$stable`/`$changed`/`$rose`/`$fell` all expand into
      // (or read from) a chain of 1-cycle latch state vars;
      // `$onehot`/`$onehot0` are a plain bit trick; `$isunknown` is a
      // constant given this encoder's pure 2-valued bitvector model.
      // Other calls (user subroutines, system tasks) are not
      // supported.
      auto & call = expr.as<CallExpression>();
      if (call.isSystemCall()
          && (call.getSubroutineName() == "$signed"
              || call.getSubroutineName() == "$unsigned")) {
        // Pure bit-pattern reinterpretation: same width, same bits,
        // only the *type* (and so, downstream, which comparison/
        // division/extension semantics apply) changes. expr.type
        // already reflects the cast's signedness for any caller that
        // inspects it (e.g. BinaryOp's op_signed check).
        auto args = call.arguments();
        if (args.empty() || !args[0]) {
          throw PonoException("SystemVerilogEncoder: "
                              + std::string(call.getSubroutineName())
                              + " with no value argument");
        }
        return expr_to_term(*args[0]);
      }
      if (call.isSystemCall() && call.getSubroutineName() == "$past") {
        auto args = call.arguments();
        if (args.empty() || !args[0]) {
          throw PonoException(
              "SystemVerilogEncoder: $past with no value argument");
        }
        Term val = expr_to_term(*args[0]);
        uint32_t n = 1;
        if (args.size() >= 2 && args[1]) {
          auto cv = args[1]->eval(eval_ctx());
          if (cv.isInteger()) {
            auto opt = cv.integer().as<uint32_t>();
            if (opt) n = *opt;
          }
        }
        if (n == 0) return val;
        return make_history_chain(val, n);
      }
      if (call.isSystemCall() && call.getSubroutineName() == "$stable") {
        auto args = call.arguments();
        if (args.empty() || !args[0]) {
          throw PonoException(
              "SystemVerilogEncoder: $stable with no value argument");
        }
        Term val = expr_to_term(*args[0]);
        Term eq = solver_->make_term(Equal, val, make_history_chain(val, 1));
        Sort bv1 = solver_->make_sort(BV, 1);
        return solver_->make_term(
            Ite, eq, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
      }
      if (call.isSystemCall() && call.getSubroutineName() == "$changed") {
        // $stable's negation: the value differs from one cycle ago,
        // over its full width (unlike $rose/$fell, which per the LRM
        // only look at bit 0).
        auto args = call.arguments();
        if (args.empty() || !args[0]) {
          throw PonoException(
              "SystemVerilogEncoder: $changed with no value argument");
        }
        Term val = expr_to_term(*args[0]);
        Term neq =
            solver_->make_term(Distinct, val, make_history_chain(val, 1));
        Sort bv1 = solver_->make_sort(BV, 1);
        return solver_->make_term(
            Ite, neq, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
      }
      if (call.isSystemCall()
          && (call.getSubroutineName() == "$rose"
              || call.getSubroutineName() == "$fell")) {
        // Per the LRM, $rose/$fell only look at bit 0 of a
        // multi-bit argument. Builds a fresh 1-cycle latch chain for
        // just that bit, via the same make_history_chain() helper
        // $past/$stable use for "value one cycle ago" -- not a chain
        // shared with any $past/$stable call on the full value.
        bool is_rose = call.getSubroutineName() == "$rose";
        auto args = call.arguments();
        if (args.empty() || !args[0]) {
          throw PonoException("SystemVerilogEncoder: "
                              + std::string(call.getSubroutineName())
                              + " with no value argument");
        }
        Term val = expr_to_term(*args[0]);
        Sort bv1 = solver_->make_sort(BV, 1);
        Term bit0 = solver_->make_term(Op(Extract, 0, 0), val);
        Term prev_bit0 = make_history_chain(bit0, 1);
        Term now_val = solver_->make_term(is_rose ? 1 : 0, bv1);
        Term prev_val = solver_->make_term(is_rose ? 0 : 1, bv1);
        Term edge =
            solver_->make_term(And,
                               solver_->make_term(Equal, bit0, now_val),
                               solver_->make_term(Equal, prev_bit0, prev_val));
        return solver_->make_term(
            Ite, edge, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
      }
      if (call.isSystemCall() && call.getSubroutineName() == "$isunknown") {
        // This encoder's SMT model is pure 2-valued bitvectors -- there
        // is no X/Z representation at all -- so nothing is ever unknown.
        return solver_->make_term(0, solver_->make_sort(BV, 1));
      }
      if (call.isSystemCall()
          && (call.getSubroutineName() == "$onehot"
              || call.getSubroutineName() == "$onehot0")) {
        // Standard power-of-two bit trick: (x & (x-1)) == 0 iff x has
        // at most one bit set; $onehot additionally requires x != 0.
        // No popcount adder needed.
        auto args = call.arguments();
        if (args.empty() || !args[0]) {
          throw PonoException("SystemVerilogEncoder: "
                              + std::string(call.getSubroutineName())
                              + " with no value argument");
        }
        Term val = expr_to_term(*args[0]);
        Term one = solver_->make_term(1, val->get_sort());
        Term minus_one = solver_->make_term(BVSub, val, one);
        Term at_most_one =
            solver_->make_term(Equal,
                               solver_->make_term(BVAnd, val, minus_one),
                               solver_->make_term(0, val->get_sort()));
        Term cond = at_most_one;
        if (call.getSubroutineName() == "$onehot") {
          Term nonzero = solver_->make_term(
              Distinct, val, solver_->make_term(0, val->get_sort()));
          cond = solver_->make_term(And, nonzero, at_most_one);
        }
        Sort bv1 = solver_->make_sort(BV, 1);
        return solver_->make_term(
            Ite, cond, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
      }
      throw PonoException("SystemVerilogEncoder: unsupported call to "
                          + std::string(call.getSubroutineName()));
    }

    default:
      throw PonoException("SystemVerilogEncoder: unsupported expression kind "
                          + to_string(static_cast<int>(expr.kind)));
  }
}

}  // namespace pono
