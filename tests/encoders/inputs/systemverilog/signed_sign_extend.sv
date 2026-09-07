// Widening a signed value must sign-extend, not zero-extend: when `a`
// (a 4-bit signed input) is negative, assigning it to the wider 8-bit
// signed `wide` must replicate the sign bit into the new upper bits,
// not fill them with zero.
module signed_sign_extend (input logic clk, input logic signed [3:0] a);

  logic signed [7:0] wide;

  assign wide = a;

  assert property (@(posedge clk) !a[3] || (wide[7:4] == 4'b1111));

endmodule
