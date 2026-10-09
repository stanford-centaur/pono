// A `disable iff` on a fairness assumption. Aborting an attempt part
// way through leaves the assumption holding over some suffixes and
// not others, and a justice condition says only "infinitely often";
// refused rather than silently promoted to unconditional.
module temporal_assume_disable_iff (
    input logic clk,
    input logic rst,
    input logic a
);

  assume property (@(posedge clk) disable iff (rst) s_eventually a);
  assert property (@(posedge clk) s_eventually a);

endmodule
