// An explicit skew samples or drives the signal somewhere inside
// the clock cycle -- a step before the edge, two time units after
// it. This encoder advances a whole cycle at a time and has no
// inside-the-cycle to place that in, so it is refused rather than
// quietly treated as the default.
module clocking_block_skew (input logic clk, input logic [3:0] d);
  logic [3:0] q;
  default clocking cb @(posedge clk);
    input #2 d;
  endclocking
  always_ff @(posedge clk) q <= cb.d;
  assert property (@(posedge clk) 1'b1);
endmodule
