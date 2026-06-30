// Tests legacy `always @*` (combinational) — encoder should detect
// no non-blocking assignments and treat the block as combinational.
// Property: y != 15.  BMC must falsify with a=12, b=3.
module legacy_always_comb (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic [3:0] y;
  always @* y = a | b;
  assert property (@(posedge clk) y != 15);
endmodule
