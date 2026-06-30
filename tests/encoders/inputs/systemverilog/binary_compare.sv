// Tests comparison operator (>=).
// Reset clears ge; one cycle later, BMC picks a, b with a >= b so ge
// becomes 1 and the property is falsified.
module binary_compare (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic ge;
  always_ff @(posedge clk) begin
    if (rst) ge <= 0;
    else ge <= a >= b;
  end
  assert property (@(posedge clk) ge != 1);
endmodule
