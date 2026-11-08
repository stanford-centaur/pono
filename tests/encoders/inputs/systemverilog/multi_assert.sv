// Tests that multiple `assert property` statements are extracted.
// All three properties are individually falsifiable.
module multi_assert (input logic clk, input logic [3:0] a);
  assert property (@(posedge clk) a != 1);
  assert property (@(posedge clk) a != 2);
  assert property (@(posedge clk) a != 3);
endmodule
