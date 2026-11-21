// Concatenation-target LHS on a plain continuous assign (`assign {hi,
// lo} = ...;`), as opposed to a concatenation-target *port connection*
// (`.port({hi, lo})`, already supported separately). Since a
// concatenation has more than one base symbol and can't be represented
// as a single LValueDesc, process_continuous_assign() special-cases it
// by splitting the RHS across each operand rather than going through
// resolve_lvalue() directly.
module concat_lhs (input logic clk,
                    input logic [3:0] a,
                    input logic [3:0] b);

  logic [3:0] hi, lo;

  assign { hi, lo } = { a, b };

  assert property (@(posedge clk) hi == a && lo == b);

endmodule
