// Tests bitwise-NOT (~).
// Reset clears x; one cycle later, free input a=5 makes x = 4'b1010 =
// 10, falsifying the property.
module unary_not (input logic clk, input logic rst, input logic [3:0] a);
  logic [3:0] x;
  always_ff @(posedge clk) begin
    if (rst) x <= 0;
    else x <= ~a;
  end
  assert property (@(posedge clk) x != 10);
endmodule
