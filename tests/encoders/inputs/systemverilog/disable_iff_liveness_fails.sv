// The paired refutation for disable_iff_liveness.sv: `c` never
// holds, so no attempt is ever aborted and the unmet eventuality is
// a genuine violation. Without this, widening the exemption to "the
// condition never holds from here on" could exempt everything and
// still look like a pass.
module disable_iff_liveness_fails (
    input logic clk
);

  logic c;
  logic b;

  initial c = 1'b0;
  initial b = 1'b0;

  always_ff @(posedge clk) c <= 1'b0;
  always_ff @(posedge clk) b <= 1'b0;

  assert property (@(posedge clk) disable iff (c) s_eventually b);

endmodule
