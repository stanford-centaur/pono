// A temporal (non-safety) `assume property`: a fairness constraint.
// `s_eventually a` at every tick is `G F a`, and a counterexample to
// any liveness property now has to keep it -- so the assertion, the
// same eventuality the assumption grants, cannot be violated by any
// fair lasso.
//
// Paired with temporal_assume_absent.sv, which is this design with
// the assumption deleted and is refuted. Without that pair, a
// fairness constraint that was quietly dropped would look the same
// as one that holds, since BMC only ever fails to find a lasso.
module temporal_assume_property (input logic clk, input logic a);

  assume property (@(posedge clk) s_eventually a);
  assert property (@(posedge clk) s_eventually a);

endmodule
