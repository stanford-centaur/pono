// Tests unary minus (-).
// Reset clears x; one cycle later, free input a=5 gives x = -5 = 11
// in 4-bit two's complement, falsifying the property.
module unary_minus (input logic clk, input logic rst, input logic [3:0] a);
  logic [3:0] x;
  always_ff @(posedge clk) begin
    if (rst) x <= 0;
    else x <= -a;
  end
  assert property (@(posedge clk) x != 11);
endmodule
