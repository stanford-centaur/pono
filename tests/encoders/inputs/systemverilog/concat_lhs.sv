// Concatenation-target LHS on a plain continuous assign (`assign {hi,
// lo} = ...;`), as opposed to a concatenation-target *port connection*
// (`.port({hi, lo})`, already supported separately). resolve_lvalue()
// has no case for ExpressionKind::Concatenation, so this silently
// failed to resolve and the whole write was dropped: hi/lo were
// declared as free/unconstrained state vars instead of being pinned.
module concat_lhs (input logic clk,
                    input logic [3:0] a,
                    input logic [3:0] b);

  logic [3:0] hi, lo;

  assign { hi, lo } = { a, b };

  assert property (@(posedge clk) hi == a && lo == b);

endmodule
