// `disable iff` on a property that stays on the LTL tableau path: an
// unbounded eventuality has nothing to re-anchor, so it cannot fall
// back to the safety encoding the way disable_iff_temporal.sv does.
//
// A liveness attempt is only ever violated by never completing, so
// its evaluation runs forever and a condition that rises *later*
// still aborts it. Here `c` is low at cycle 0 and high from cycle 1
// on, so every attempt is aborted even though `b` never occurs, and
// nothing is violated. Exempting only the attempt's own start cycle
// reports a violation instead -- the assertion made stronger than
// written. disable_iff_liveness_fails.sv is the same design with `c`
// held low, where the violation is real.
module disable_iff_liveness (
    input logic clk
);

  logic c;
  logic b;

  initial c = 1'b0;
  initial b = 1'b0;

  always_ff @(posedge clk) c <= 1'b1;
  always_ff @(posedge clk) b <= 1'b0;

  assert property (@(posedge clk) disable iff (c) s_eventually b);

endmodule
