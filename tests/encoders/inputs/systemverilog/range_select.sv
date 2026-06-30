// Tests a constant range-select expression a[7:4].
// Reset clears x; one cycle later, the upper nibble of the free
// input a can be 10, falsifying the property.
module range_select (input logic clk, input logic rst, input logic [7:0] a);
  logic [3:0] x;
  always_ff @(posedge clk) begin
    if (rst) x <= 0;
    else x <= a[7:4];
  end
  assert property (@(posedge clk) x != 10);
endmodule
