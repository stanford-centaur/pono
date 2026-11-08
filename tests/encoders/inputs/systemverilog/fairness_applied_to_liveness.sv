// The other half of fairness_not_applied_to_safety.sv: the same
// design, asked as a liveness property, where the fairness
// assumption does reach. Reaching `dead` would mean `a` never rises
// again, so this fails exactly when a fair trace can get there --
// and it is proved, so none can.
//
// The pair is the point: one design, one assumption, two verdicts
// that disagree about whether `dead` happens, and both correct for
// the class of trace each is asking about.
module fairness_applied_to_liveness (
    input logic clk,
    input logic a,
    input logic go
);

  logic dead;
  initial dead = 1'b0;
  always_ff @(posedge clk) if (go) dead <= 1'b1;

  // Once dead, the fairness condition can never be met again.
  assume property (@(posedge clk) dead |-> !a);
  assume property (@(posedge clk) s_eventually a);

  assert property (@(posedge clk) always (dead |-> s_eventually a));

endmodule
