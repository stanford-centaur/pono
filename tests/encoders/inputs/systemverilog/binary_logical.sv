// Tests logical operators (&&, ||).
// Reset clears y; one cycle later, free a=1 && b=1 sets y to 1 and
// falsifies the property.
module binary_logical (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b
);
  logic y;
  always_ff @(posedge clk) begin
    if (rst) y <= 0;
    else y <= a && b;
  end
  assert property (@(posedge clk) y != 1);
endmodule
