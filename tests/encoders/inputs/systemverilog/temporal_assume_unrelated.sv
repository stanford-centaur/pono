// The other half of the pair, and the one that catches an
// assumption that constrains too much rather than too little: the
// assumption names `a` and the assertion an unrelated free `b`, so
// it must still be refuted. An assumption whose tableau were
// unsatisfiable would empty the model and "prove" everything,
// which temporal_assume_property.sv alone cannot tell apart from
// the fairness constraint working.
module temporal_assume_unrelated (
    input logic clk,
    input logic a,
    input logic b
);

  assume property (@(posedge clk) s_eventually a);
  assert property (@(posedge clk) s_eventually b);

endmodule
