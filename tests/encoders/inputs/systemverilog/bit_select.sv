// Tests a constant bit-select expression a[3].
// Property: x != 1.  BMC must falsify when a free input has bit 3 set.
module bit_select (input logic clk, input logic [7:0] a);
  logic x;
  always_ff @(posedge clk) x <= a[3];
  assert property (@(posedge clk) x != 1);
endmodule
