// A `cover property` whose goal is a genuine eventuality. A
// sequence-shaped goal works -- a bare sequence is a per-cycle check
// anchored at its own attempt, not an eventuality -- but one built
// from `s_eventually`, an unbounded delay or a counted repetition
// has to run through the LTL tableau, where reachability duality
// needs the verdict reported inverted and Result carries no flag
// saying so.
//
// Intended: the goal is reachable (`a` is free), so the dual assert
// `G !a` has a fair lasso and the cover succeeds.
module cover_temporal_goal (input logic clk, input logic a);

  cover property (@(posedge clk) s_eventually a);

endmodule
