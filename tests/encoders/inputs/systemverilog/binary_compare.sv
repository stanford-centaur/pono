// Tests comparison operator (>=).
// Property: ge != 1.  BMC must falsify within 2 steps because the
// register `ge` gets a >= b, and a, b are free inputs.
module binary_compare (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic ge;
  always_ff @(posedge clk) ge <= a >= b;
  assert property (@(posedge clk) ge != 1);
endmodule
