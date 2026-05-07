// Tests bit replication {{N{x}}}.
// Property: x != 255.  BMC must falsify with a free 2-bit input a = 3
// (since {{4{2'b11}}} = 8'b11111111 = 255).
module replication (input logic clk, input logic [1:0] a);
  logic [7:0] x;
  always_ff @(posedge clk) x <= {4{a}};
  assert property (@(posedge clk) x != 255);
endmodule
