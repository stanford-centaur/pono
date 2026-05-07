// Tests unary minus (-).
// Property: x != 11.  BMC must falsify when a == 5 because in 4-bit
// two's complement -5 = 11.
module unary_minus (input logic clk, input logic [3:0] a);
  logic [3:0] x;
  always_ff @(posedge clk) x <= -a;
  assert property (@(posedge clk) x != 11);
endmodule
