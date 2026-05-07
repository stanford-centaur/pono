// Tests bitwise-NOT (~).
// Property: x != 10.  BMC must falsify when a == 5
// (since ~4'd5 = 4'b1010 = 10).
module unary_not (input logic clk, input logic [3:0] a);
  logic [3:0] x;
  always_ff @(posedge clk) x <= ~a;
  assert property (@(posedge clk) x != 10);
endmodule
