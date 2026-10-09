// Relational operators on signed operands must compare as two's-
// complement signed values, not as raw unsigned bit patterns: when
// a == -1 (4'b1111) and b == 1, `a < b` must hold under signed
// semantics even though the raw bit pattern of a (15) is numerically
// larger than b's (1) under unsigned semantics.
module signed_compare (input logic clk,
                       input logic signed [3:0] a,
                       input logic signed [3:0] b);

  assert property (@(posedge clk) !(a == -4'sd1 && b == 4'sd1) || (a < b));

endmodule
