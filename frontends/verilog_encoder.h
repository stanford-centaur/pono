/*********************                                                        */
/*! \file verilog_encoder.h
** \verbatim
** Top contributors (to current version):
**   Codex
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Limited Verilog 2005 frontend (handwritten parser).
**
**/

#pragma once

#include <cctype>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "core/ts.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"

namespace pono {

class VerilogEncoder
{
 public:
  VerilogEncoder(const std::string & filename, TransitionSystem & ts)
      : ts_(ts), solver_(ts.solver())
  {
    std::ifstream in(filename);
    if (!in) {
      throw PonoException("Could not open Verilog file: " + filename);
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    source_ = ss.str();
    lex();
    parse();
  }

 private:
  enum class TokenKind
  {
    End,
    Identifier,
    Number,
    Symbol
  };

  struct Token
  {
    TokenKind kind;
    std::string text;
    size_t line;
    size_t col;
  };

  struct Expr
  {
    smt::Term term;
    uint64_t width;
    bool is_bool;
  };

  enum class SignalKind
  {
    Input,
    Reg,
    Wire,
    Output
  };

  struct SignalInfo
  {
    SignalKind kind;
    uint64_t width;
    smt::Term term;  // valid for inputs/regs
    Expr expr;       // valid for wires/outputs with assign
    bool has_expr = false;
  };

  struct Statement
  {
    enum class Kind
    {
      Block,
      Assign,
      If
    };

    Kind kind;
    std::vector<Statement> children;  // block statements
    std::string lhs;                  // assign target (identifier)
    Expr rhs;                         // assign rhs
    Expr cond;                        // if condition
    std::unique_ptr<Statement> then_s;
    std::unique_ptr<Statement> else_s;
  };

  TransitionSystem & ts_;
  const smt::SmtSolver & solver_;
  std::string source_;
  std::vector<Token> tokens_;
  size_t pos_ = 0;

  std::unordered_map<std::string, SignalInfo> signals_;
  std::unordered_map<std::string, Expr> wire_exprs_;
  std::unordered_map<std::string, smt::Term> state_terms_;
  std::unordered_set<std::string> assigned_regs_;

  void lex()
  {
    size_t i = 0;
    size_t line = 1;
    size_t col = 1;
    auto push = [&](TokenKind kind, const std::string & text) {
      tokens_.push_back({ kind, text, line, col });
    };
    while (i < source_.size()) {
      char c = source_[i];
      if (c == '\n') {
        ++line;
        col = 1;
        ++i;
        continue;
      }
      if (std::isspace(static_cast<unsigned char>(c))) {
        ++col;
        ++i;
        continue;
      }
      if (c == '/' && i + 1 < source_.size() && source_[i + 1] == '/') {
        i += 2;
        col += 2;
        while (i < source_.size() && source_[i] != '\n') {
          ++i;
          ++col;
        }
        continue;
      }
      if (c == '/' && i + 1 < source_.size() && source_[i + 1] == '*') {
        i += 2;
        col += 2;
        while (i + 1 < source_.size()
               && !(source_[i] == '*' && source_[i + 1] == '/')) {
          if (source_[i] == '\n') {
            ++line;
            col = 1;
          } else {
            ++col;
          }
          ++i;
        }
        if (i + 1 < source_.size()) {
          i += 2;
          col += 2;
        }
        continue;
      }
      if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
        size_t start = i;
        size_t start_col = col;
        while (i < source_.size()
               && (std::isalnum(static_cast<unsigned char>(source_[i]))
                   || source_[i] == '_')) {
          ++i;
          ++col;
        }
        tokens_.push_back(
            { TokenKind::Identifier, source_.substr(start, i - start), line,
              start_col });
        continue;
      }
      if (std::isdigit(static_cast<unsigned char>(c))) {
        size_t start = i;
        size_t start_col = col;
        while (i < source_.size()
               && (std::isalnum(static_cast<unsigned char>(source_[i]))
                   || source_[i] == '\'')) {
          ++i;
          ++col;
        }
        tokens_.push_back(
            { TokenKind::Number, source_.substr(start, i - start), line,
              start_col });
        continue;
      }
      // operators and punctuation
      auto two = (i + 1 < source_.size()) ? source_.substr(i, 2) : "";
      auto three = (i + 2 < source_.size()) ? source_.substr(i, 3) : "";
      if (three == "<<<" || three == ">>>" || three == "&&&") {
        push(TokenKind::Symbol, three);
        i += 3;
        col += 3;
        continue;
      }
      if (two == "==" || two == "!=" || two == "<=" || two == ">="
          || two == "&&" || two == "||" || two == "<<" || two == ">>"
          || two == "<=" || two == ">=" || two == "<>") {
        push(TokenKind::Symbol, two);
        i += 2;
        col += 2;
        continue;
      }
      push(TokenKind::Symbol, std::string(1, c));
      ++i;
      ++col;
    }
    tokens_.push_back({ TokenKind::End, "", line, col });
  }

  const Token & peek() const { return tokens_.at(pos_); }

  const Token & advance()
  {
    if (pos_ < tokens_.size()) {
      ++pos_;
    }
    return tokens_.at(pos_ - 1);
  }

  bool match_symbol(const std::string & sym)
  {
    if (peek().kind == TokenKind::Symbol && peek().text == sym) {
      advance();
      return true;
    }
    return false;
  }

  bool match_keyword(const std::string & kw)
  {
    if (peek().kind == TokenKind::Identifier && peek().text == kw) {
      advance();
      return true;
    }
    return false;
  }

  void expect_symbol(const std::string & sym)
  {
    if (!match_symbol(sym)) {
      throw parse_error("Expected symbol '" + sym + "'");
    }
  }

  std::string expect_identifier(const std::string & context)
  {
    if (peek().kind != TokenKind::Identifier) {
      throw parse_error("Expected identifier in " + context);
    }
    return advance().text;
  }

  PonoException parse_error(const std::string & msg) const
  {
    std::ostringstream ss;
    ss << "Verilog parse error at " << peek().line << ":" << peek().col
       << ": " << msg;
    return PonoException(ss.str());
  }

  void parse()
  {
    if (!match_keyword("module")) {
      throw parse_error("Expected 'module'");
    }
    expect_identifier("module name");
    // consume optional port list
    if (match_symbol("(")) {
      while (!match_symbol(")")) {
        if (peek().kind == TokenKind::End) {
          throw parse_error("Unterminated port list");
        }
        advance();
      }
    }
    expect_symbol(";");

    while (!match_keyword("endmodule")) {
      if (peek().kind == TokenKind::End) {
        throw parse_error("Missing 'endmodule'");
      }
      if (match_keyword("input")) {
        parse_decl(SignalKind::Input);
      } else if (match_keyword("output")) {
        parse_output_decl();
      } else if (match_keyword("reg")) {
        parse_decl(SignalKind::Reg);
      } else if (match_keyword("wire")) {
        parse_decl(SignalKind::Wire);
      } else if (match_keyword("assign")) {
        parse_assign();
      } else if (match_keyword("always")) {
        parse_always();
      } else {
        throw parse_error("Unsupported construct");
      }
    }
  }

  void parse_decl(SignalKind kind)
  {
    uint64_t width = 1;
    if (match_symbol("[")) {
      auto msb = parse_int_constant();
      expect_symbol(":");
      auto lsb = parse_int_constant();
      expect_symbol("]");
      width = width_from_range(msb, lsb);
    }
    while (true) {
      std::string name = expect_identifier("declaration");
      declare_signal(name, kind, width);
      if (!match_symbol(",")) {
        break;
      }
    }
    expect_symbol(";");
  }

  void parse_output_decl()
  {
    SignalKind kind = SignalKind::Wire;
    if (match_keyword("reg")) {
      kind = SignalKind::Reg;
    } else if (match_keyword("wire")) {
      kind = SignalKind::Wire;
    }
    uint64_t width = 1;
    if (match_symbol("[")) {
      auto msb = parse_int_constant();
      expect_symbol(":");
      auto lsb = parse_int_constant();
      expect_symbol("]");
      width = width_from_range(msb, lsb);
    }
    while (true) {
      std::string name = expect_identifier("output declaration");
      declare_signal(name, kind, width);
      if (!match_symbol(",")) {
        break;
      }
    }
    expect_symbol(";");
  }

  void declare_signal(const std::string & name,
                      SignalKind kind,
                      uint64_t width)
  {
    if (signals_.count(name)) {
      throw PonoException("Duplicate signal declaration: " + name);
    }
    SignalInfo info;
    info.kind = kind;
    info.width = width;
    if (kind == SignalKind::Input) {
      info.term = ts_.make_inputvar(name, solver_->make_sort(smt::BV, width));
    } else if (kind == SignalKind::Reg) {
      info.term = ts_.make_statevar(name, solver_->make_sort(smt::BV, width));
      state_terms_[name] = info.term;
    }
    signals_.emplace(name, info);
  }

  void parse_assign()
  {
    std::string lhs = expect_identifier("assign");
    expect_symbol("=");
    Expr rhs = parse_expr();
    expect_symbol(";");
    auto it = signals_.find(lhs);
    if (it == signals_.end()) {
      throw PonoException("Assign to undeclared signal: " + lhs);
    }
    if (it->second.kind == SignalKind::Reg) {
      throw PonoException("Continuous assign to reg not supported: " + lhs);
    }
    it->second.expr = ensure_width(rhs, it->second.width);
    it->second.has_expr = true;
    wire_exprs_[lhs] = it->second.expr;
  }

  void parse_always()
  {
    if (!match_symbol("@")) {
      throw parse_error("Expected '@' after always");
    }
    expect_symbol("(");
    bool saw_edge = false;
    while (!match_symbol(")")) {
      if (peek().kind == TokenKind::End) {
        throw parse_error("Unterminated sensitivity list");
      }
      if (match_keyword("posedge") || match_keyword("negedge")) {
        saw_edge = true;
        expect_identifier("clock signal");
      } else if (match_keyword("or")) {
        continue;
      } else {
        advance();
      }
    }
    if (!saw_edge) {
      throw PonoException("Only edge-triggered always blocks are supported");
    }
    Statement stmt = parse_statement();
    build_next_from_statement(stmt);
  }

  Statement parse_statement()
  {
    if (match_keyword("begin")) {
      Statement st;
      st.kind = Statement::Kind::Block;
      while (!match_keyword("end")) {
        if (peek().kind == TokenKind::End) {
          throw parse_error("Unterminated begin/end block");
        }
        st.children.push_back(parse_statement());
      }
      return st;
    }
    if (match_keyword("if")) {
      expect_symbol("(");
      Expr cond = parse_expr();
      expect_symbol(")");
      Statement st;
      st.kind = Statement::Kind::If;
      st.cond = to_bool(cond);
      st.then_s.reset(new Statement(parse_statement()));
      if (match_keyword("else")) {
        st.else_s.reset(new Statement(parse_statement()));
      }
      return st;
    }
    // assignment
    std::string lhs = expect_identifier("assignment");
    if (!match_symbol("<=")) {
      throw parse_error("Only nonblocking assignments are supported");
    }
    Expr rhs = parse_expr();
    expect_symbol(";");
    Statement st;
    st.kind = Statement::Kind::Assign;
    st.lhs = lhs;
    st.rhs = rhs;
    return st;
  }

  void build_next_from_statement(const Statement & stmt)
  {
    smt::Term cond_true = solver_->make_term(true);
    std::unordered_map<std::string, Expr> env;
    for (const auto & it : state_terms_) {
      env[it.first] = { it.second, signals_.at(it.first).width, false };
    }
    apply_statement(stmt, cond_true, env);
    for (const auto & it : env) {
      const std::string & name = it.first;
      const Expr & next_expr = it.second;
      if (!assigned_regs_.count(name)) {
        ts_.assign_next(state_terms_.at(name),
                        ensure_width(next_expr, signals_.at(name).width).term);
        assigned_regs_.insert(name);
      } else {
        throw PonoException("Multiple always blocks assign to reg: " + name);
      }
    }
  }

  void apply_statement(const Statement & stmt,
                       const smt::Term & cond,
                       std::unordered_map<std::string, Expr> & env)
  {
    if (stmt.kind == Statement::Kind::Block) {
      for (const auto & child : stmt.children) {
        apply_statement(child, cond, env);
      }
      return;
    }
    if (stmt.kind == Statement::Kind::If) {
      Expr c = to_bool(stmt.cond);
      smt::Term then_c =
          solver_->make_term(smt::And, cond, c.term);
      smt::Term else_c =
          solver_->make_term(smt::And, cond,
                             solver_->make_term(smt::Not, c.term));
      if (stmt.then_s) {
        apply_statement(*stmt.then_s, then_c, env);
      }
      if (stmt.else_s) {
        apply_statement(*stmt.else_s, else_c, env);
      }
      return;
    }
    if (stmt.kind == Statement::Kind::Assign) {
      auto it = signals_.find(stmt.lhs);
      if (it == signals_.end() || it->second.kind != SignalKind::Reg) {
        throw PonoException("Assignment target is not a reg: " + stmt.lhs);
      }
      Expr rhs = ensure_width(stmt.rhs, it->second.width);
      Expr curr = env.at(stmt.lhs);
      smt::Term ite =
          solver_->make_term(smt::Ite, cond, rhs.term, curr.term);
      env[stmt.lhs] = { ite, it->second.width, false };
      return;
    }
    throw PonoException("Unhandled statement kind");
  }

  uint64_t parse_int_constant()
  {
    if (peek().kind != TokenKind::Number) {
      throw parse_error("Expected integer constant");
    }
    std::string text = advance().text;
    auto tick = text.find('\'');
    if (tick == std::string::npos) {
      return std::stoull(text);
    }
    std::string size = text.substr(0, tick);
    std::string base_and_val = text.substr(tick + 1);
    if (base_and_val.empty()) {
      throw parse_error("Malformed sized literal");
    }
    char base = base_and_val[0];
    std::string val = base_and_val.substr(1);
    for (char ch : val) {
      if (ch == 'x' || ch == 'z' || ch == 'X' || ch == 'Z') {
        throw PonoException("x/z literals are not supported");
      }
    }
    uint64_t v = 0;
    int radix = 10;
    if (base == 'b' || base == 'B') {
      radix = 2;
    } else if (base == 'h' || base == 'H') {
      radix = 16;
    } else if (base == 'o' || base == 'O') {
      radix = 8;
    } else if (base == 'd' || base == 'D') {
      radix = 10;
    } else {
      throw parse_error("Unknown literal base");
    }
    v = std::stoull(val, nullptr, radix);
    return v;
  }

  std::pair<uint64_t, uint64_t> parse_number_literal()
  {
    if (peek().kind != TokenKind::Number) {
      throw parse_error("Expected number literal");
    }
    std::string text = advance().text;
    auto tick = text.find('\'');
    if (tick == std::string::npos) {
      return { std::stoull(text), 32 };
    }
    uint64_t width = std::stoull(text.substr(0, tick));
    std::string base_and_val = text.substr(tick + 1);
    if (base_and_val.empty()) {
      throw parse_error("Malformed sized literal");
    }
    char base = base_and_val[0];
    std::string val = base_and_val.substr(1);
    for (char ch : val) {
      if (ch == 'x' || ch == 'z' || ch == 'X' || ch == 'Z') {
        throw PonoException("x/z literals are not supported");
      }
    }
    int radix = 10;
    if (base == 'b' || base == 'B') {
      radix = 2;
    } else if (base == 'h' || base == 'H') {
      radix = 16;
    } else if (base == 'o' || base == 'O') {
      radix = 8;
    } else if (base == 'd' || base == 'D') {
      radix = 10;
    } else {
      throw parse_error("Unknown literal base");
    }
    uint64_t value = std::stoull(val, nullptr, radix);
    return { value, width };
  }

  uint64_t width_from_range(uint64_t msb, uint64_t lsb)
  {
    if (msb >= lsb) {
      return msb - lsb + 1;
    }
    return lsb - msb + 1;
  }

  Expr parse_expr() { return parse_logical_or(); }

  Expr parse_logical_or()
  {
    Expr left = parse_logical_and();
    while (match_symbol("||")) {
      Expr right = parse_logical_and();
      left = { solver_->make_term(smt::Or, to_bool(left).term,
                                  to_bool(right).term),
               1, true };
    }
    return left;
  }

  Expr parse_logical_and()
  {
    Expr left = parse_bitor();
    while (match_symbol("&&")) {
      Expr right = parse_bitor();
      left = { solver_->make_term(smt::And, to_bool(left).term,
                                  to_bool(right).term),
               1, true };
    }
    return left;
  }

  Expr parse_bitor()
  {
    Expr left = parse_bitxor();
    while (match_symbol("|")) {
      Expr right = parse_bitxor();
      left = bitwise_binop(smt::BVOr, left, right);
    }
    return left;
  }

  Expr parse_bitxor()
  {
    Expr left = parse_bitand();
    while (match_symbol("^")) {
      Expr right = parse_bitand();
      left = bitwise_binop(smt::BVXor, left, right);
    }
    return left;
  }

  Expr parse_bitand()
  {
    Expr left = parse_equality();
    while (match_symbol("&")) {
      Expr right = parse_equality();
      left = bitwise_binop(smt::BVAnd, left, right);
    }
    return left;
  }

  Expr parse_equality()
  {
    Expr left = parse_relational();
    while (true) {
      if (match_symbol("==")) {
        Expr right = parse_relational();
        if (left.is_bool || right.is_bool) {
          left = { solver_->make_term(smt::Equal, to_bool(left).term,
                                      to_bool(right).term),
                   1, true };
        } else {
          auto l = ensure_width(left, std::max(left.width, right.width));
          auto r = ensure_width(right, l.width);
          left = { solver_->make_term(smt::Equal, l.term, r.term), 1, true };
        }
      } else if (match_symbol("!=")) {
        Expr right = parse_relational();
        smt::Term eq;
        if (left.is_bool || right.is_bool) {
          eq = solver_->make_term(smt::Equal, to_bool(left).term,
                                  to_bool(right).term);
        } else {
          auto l = ensure_width(left, std::max(left.width, right.width));
          auto r = ensure_width(right, l.width);
          eq = solver_->make_term(smt::Equal, l.term, r.term);
        }
        left = { solver_->make_term(smt::Not, eq), 1, true };
      } else {
        break;
      }
    }
    return left;
  }

  Expr parse_relational()
  {
    Expr left = parse_shift();
    while (true) {
      if (match_symbol("<")) {
        Expr right = parse_shift();
        auto l = ensure_width(left, std::max(left.width, right.width));
        auto r = ensure_width(right, l.width);
        left = { solver_->make_term(smt::BVUlt, l.term, r.term), 1, true };
      } else if (match_symbol("<=")) {
        Expr right = parse_shift();
        auto l = ensure_width(left, std::max(left.width, right.width));
        auto r = ensure_width(right, l.width);
        left = { solver_->make_term(smt::BVUle, l.term, r.term), 1, true };
      } else if (match_symbol(">")) {
        Expr right = parse_shift();
        auto l = ensure_width(left, std::max(left.width, right.width));
        auto r = ensure_width(right, l.width);
        left = { solver_->make_term(smt::BVUgt, l.term, r.term), 1, true };
      } else if (match_symbol(">=")) {
        Expr right = parse_shift();
        auto l = ensure_width(left, std::max(left.width, right.width));
        auto r = ensure_width(right, l.width);
        left = { solver_->make_term(smt::BVUge, l.term, r.term), 1, true };
      } else {
        break;
      }
    }
    return left;
  }

  Expr parse_shift()
  {
    Expr left = parse_add();
    while (true) {
      if (match_symbol("<<")) {
        Expr right = parse_add();
        Expr r = ensure_width(right, left.width);
        left = { solver_->make_term(smt::BVShl, left.term, r.term), left.width,
                 false };
      } else if (match_symbol(">>")) {
        Expr right = parse_add();
        Expr r = ensure_width(right, left.width);
        left = { solver_->make_term(smt::BVShrl, left.term, r.term),
                 left.width, false };
      } else {
        break;
      }
    }
    return left;
  }

  Expr parse_add()
  {
    Expr left = parse_mul();
    while (true) {
      if (match_symbol("+")) {
        Expr right = parse_mul();
        left = arith_binop(smt::BVAdd, left, right);
      } else if (match_symbol("-")) {
        Expr right = parse_mul();
        left = arith_binop(smt::BVSub, left, right);
      } else {
        break;
      }
    }
    return left;
  }

  Expr parse_mul()
  {
    Expr left = parse_unary();
    while (match_symbol("*")) {
      Expr right = parse_unary();
      left = arith_binop(smt::BVMul, left, right);
    }
    return left;
  }

  Expr parse_unary()
  {
    if (match_symbol("!")) {
      Expr e = parse_unary();
      return { solver_->make_term(smt::Not, to_bool(e).term), 1, true };
    }
    if (match_symbol("~")) {
      Expr e = parse_unary();
      return { solver_->make_term(smt::BVNot, ensure_bv(e).term), e.width,
               false };
    }
    if (match_symbol("-")) {
      Expr e = parse_unary();
      return { solver_->make_term(smt::BVNeg, ensure_bv(e).term), e.width,
               false };
    }
    return parse_primary();
  }

  Expr parse_primary()
  {
    if (match_symbol("(")) {
      Expr e = parse_expr();
      expect_symbol(")");
      return e;
    }
    if (peek().kind == TokenKind::Number) {
      auto lit = parse_number_literal();
      uint64_t val = lit.first;
      uint64_t width = lit.second ? lit.second : 32;
      smt::Term t = solver_->make_term(val, solver_->make_sort(smt::BV, width));
      return { t, width, false };
    }
    std::string id = expect_identifier("expression");
    auto it = signals_.find(id);
    if (it == signals_.end()) {
      throw PonoException("Unknown identifier: " + id);
    }
    Expr base = resolve_signal_expr(id);
    if (match_symbol("[")) {
      uint64_t msb = parse_int_constant();
      if (match_symbol(":")) {
        uint64_t lsb = parse_int_constant();
        expect_symbol("]");
        uint64_t hi = std::max(msb, lsb);
        uint64_t lo = std::min(msb, lsb);
        uint64_t width = width_from_range(msb, lsb);
        smt::Term t = solver_->make_term(
            smt::Op(smt::PrimOp::Extract, hi, lo), base.term);
        return { t, width, false };
      }
      expect_symbol("]");
      smt::Term t =
          solver_->make_term(smt::Op(smt::PrimOp::Extract, msb, msb),
                             base.term);
      return { t, 1, false };
    }
    return base;
  }

  Expr resolve_signal_expr(const std::string & name)
  {
    auto it = signals_.find(name);
    if (it == signals_.end()) {
      throw PonoException("Unknown signal: " + name);
    }
    if (it->second.kind == SignalKind::Input
        || it->second.kind == SignalKind::Reg) {
      return { it->second.term, it->second.width, false };
    }
    if (it->second.has_expr) {
      return it->second.expr;
    }
    throw PonoException("Wire/output used before assignment: " + name);
  }

  Expr ensure_width(const Expr & e, uint64_t width)
  {
    Expr bv = ensure_bv(e);
    if (bv.width == width) {
      return bv;
    }
    if (bv.width < width) {
    smt::Term ext = solver_->make_term(
        smt::Op(smt::PrimOp::Zero_Extend, width - bv.width), bv.term);
      return { ext, width, false };
    }
    smt::Term trunc =
        solver_->make_term(smt::Op(smt::PrimOp::Extract, width - 1, 0),
                           bv.term);
    return { trunc, width, false };
  }

  Expr ensure_bv(const Expr & e)
  {
    if (!e.is_bool) {
      return e;
    }
    smt::Term one =
        solver_->make_term(1, solver_->make_sort(smt::BV, 1));
    smt::Term zero =
        solver_->make_term(0, solver_->make_sort(smt::BV, 1));
    smt::Term t = solver_->make_term(smt::Ite, e.term, one, zero);
    return { t, 1, false };
  }

  Expr to_bool(const Expr & e)
  {
    if (e.is_bool) {
      return e;
    }
    smt::Term zero =
        solver_->make_term(0, solver_->make_sort(smt::BV, e.width));
    smt::Term eq = solver_->make_term(smt::Equal, e.term, zero);
    return { solver_->make_term(smt::Not, eq), 1, true };
  }

  Expr bitwise_binop(smt::PrimOp op, Expr left, Expr right)
  {
    uint64_t width = std::max(left.width, right.width);
    left = ensure_width(left, width);
    right = ensure_width(right, width);
    return { solver_->make_term(op, left.term, right.term), width, false };
  }

  Expr arith_binop(smt::PrimOp op, Expr left, Expr right)
  {
    uint64_t width = std::max(left.width, right.width);
    left = ensure_width(left, width);
    right = ensure_width(right, width);
    return { solver_->make_term(op, left.term, right.term), width, false };
  }
};

}  // namespace pono
