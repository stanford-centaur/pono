// Tests bitwise operators (&, |, ^).
// Reset clears x; one cycle later, free inputs with a | b == 15
// (e.g. a=12, b=3) falsify the property.
module binary_bitwise (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic [3:0] x;
  always_ff @(posedge clk) begin
    if (rst) x <= 0;
    else x <= a | b;
  end
  assert property (@(posedge clk) x != 15);
endmodule
