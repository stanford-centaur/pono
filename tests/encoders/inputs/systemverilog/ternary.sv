// Tests the ?: ternary operator.
// Reset clears x; one cycle later, sel=1 selects 9 and falsifies the
// property.
module ternary (input logic clk, input logic rst, input logic sel);
  logic [3:0] x;
  always_ff @(posedge clk) begin
    if (rst) x <= 0;
    else x <= sel ? 9 : 2;
  end
  assert property (@(posedge clk) x != 9);
endmodule
