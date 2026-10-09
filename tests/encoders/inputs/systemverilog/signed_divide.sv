// Division on signed operands must truncate toward zero using signed
// semantics, not treat the operands' raw bit patterns as unsigned:
// -8 / 3 == -2 (truncating toward zero) under signed division, but
// 248 / 3 == 82 under unsigned division of the same bit patterns.
module signed_divide (input logic clk,
                      input logic signed [7:0] a,
                      input logic signed [7:0] b);

  assert property (@(posedge clk) !(a == -8'sd8 && b == 8'sd3)
                    || (a / b == -8'sd2));

endmodule
