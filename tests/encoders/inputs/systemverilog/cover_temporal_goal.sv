// A `cover property` whose goal is a genuine eventuality, built from
// `s_eventually` rather than from a bare sequence. A sequence-shaped
// goal is a per-cycle check anchored at its own attempt; this one has
// to run through the LTL tableau, where reachability duality reports
// the verdict inverted.
//
// `a` is free, so the goal is reachable and the cover succeeds.
module cover_temporal_goal (input logic clk, input logic a);

  cover property (@(posedge clk) s_eventually a);

endmodule
