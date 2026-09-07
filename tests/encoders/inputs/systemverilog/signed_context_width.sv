// The result of `a + b` is computed at the operands' own natural
// (4-bit) width first, then widened to fit the wider context
// (`wide_sum`'s 16-bit assignment target) -- that final widening must
// sign-extend, not zero-extend, or a negative sum becomes a large
// positive one purely from being assigned into a wider variable.
module signed_context_width (input logic clk,
                             input logic signed [3:0] a,
                             input logic signed [3:0] b);

  logic signed [15:0] wide_sum;

  assign wide_sum = a + b;

  assert property (@(posedge clk) !(a == -4'sd1 && b == -4'sd1)
                    || (wide_sum == -16'sd2));

endmodule
