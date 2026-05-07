// Tests left-shift operator (<<).
// Property: x != 12.  BMC must falsify when free input a == 3,
// since 4'd3 << 2 == 4'b1100 == 12.
module shift (input logic clk, input logic [3:0] a);
  logic [3:0] x;
  always_ff @(posedge clk) x <= a << 2;
  assert property (@(posedge clk) x != 12);
endmodule
