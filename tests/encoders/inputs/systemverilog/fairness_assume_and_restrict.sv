// `assume` and `restrict` in one design. Only simulation tells the
// two apart, so both land in the same pool of fairness conditions.
// The assertion is `s_eventually b` written so that only the
// `restrict` can grant it -- `a && !a` is false at every cycle, so
// the `assume` on `a` contributes nothing to it and the property
// proves only if the restrict's condition reached it too.
module fairness_assume_and_restrict (
    input logic clk,
    input logic a,
    input logic b
);

  assume property (@(posedge clk) s_eventually a);
  restrict property (@(posedge clk) s_eventually b);
  assert property (@(posedge clk) s_eventually ((a && !a) || b));

endmodule
