// Tests a constant range-select expression a[7:4].
// Property: x != 10.  BMC must falsify with a free input whose upper
// nibble is 10 (1010).
module range_select (input logic clk, input logic [7:0] a);
  logic [3:0] x;
  always_ff @(posedge clk) x <= a[7:4];
  assert property (@(posedge clk) x != 10);
endmodule
