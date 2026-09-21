// A fairness assumption alongside a safety property. The assumption
// rules out infinite traces only, and a safety counterexample is a
// finite prefix, so there is no justice condition that can exclude
// one -- telling a prefix with a fair continuation from one without
// needs a fair-reachability check this encoder has no notion of.
// Refused rather than reporting a counterexample the assumption may
// have excluded.
module temporal_assume_with_safety (input logic clk, input logic a);

  assume property (@(posedge clk) s_eventually a);
  assert property (@(posedge clk) 1'b1);

endmodule
