// A temporal (non-safety) `assume property` -- `s_eventually a`
// needs its own fairness-constraint machinery (assuming a GF
// condition rather than proving one), which this encoder doesn't
// build. Dropping it silently would leave the model less constrained
// than the source describes, risking a spurious counterexample.
module temporal_assume_property (input logic clk, input logic a);
  assume property (@(posedge clk) s_eventually a);
  assert property (@(posedge clk) 1'b1);
endmodule
