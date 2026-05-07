// Tests logical operators (&&, ||).
// Property: y != 1.  BMC must falsify when (a && b) is true OR a is
// nonzero OR-ing produces nonzero result.
module binary_logical (
    input logic clk,
    input logic a,
    input logic b
);
  logic y;
  always_ff @(posedge clk) y <= a && b;
  assert property (@(posedge clk) y != 1);
endmodule
