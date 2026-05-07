// Tests `always_comb` blocks for combinational logic.
// Property: y != 15.  BMC must falsify with a=12, b=3.
module always_comb_test (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic [3:0] y;
  always_comb y = a | b;
  assert property (@(posedge clk) y != 15);
endmodule
