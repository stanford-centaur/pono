// Tests a constant bit-select expression a[3].
// Reset clears x; one cycle later, when bit 3 of the free input a is
// 1, x reflects that and the property is falsified.
module bit_select (input logic clk, input logic rst, input logic [7:0] a);
  logic x;
  always_ff @(posedge clk) begin
    if (rst) x <= 0;
    else x <= a[3];
  end
  assert property (@(posedge clk) x != 1);
endmodule
