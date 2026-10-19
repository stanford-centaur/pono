// A `disable iff` on a fairness assumption. The assumption holds over
// the suffixes the condition does not abort, which a justice
// condition -- saying only "infinitely often" -- cannot express on
// its own.
//
// `rst` is tied low, so nothing is ever aborted and the assumption is
// the plain fairness constraint that grants the assertion.
module temporal_assume_disable_iff (
    input logic clk,
    input logic a
);

  logic rst;
  assign rst = 1'b0;

  assume property (@(posedge clk) disable iff (rst) s_eventually a);
  assert property (@(posedge clk) s_eventually a);

endmodule
