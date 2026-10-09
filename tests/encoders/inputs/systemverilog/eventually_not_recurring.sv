// Distinguishes "eventually" from "infinitely often", which is what
// the implicit per-cycle closure on a property expression amounts to.
//
// `tick` saturates at 15, so (tick == 3) is true at exactly one cycle
// of the whole trace and never again.  Per the LRM a property
// expression is evaluated at every clock tick, so this assertion means
// `always (s_eventually (tick == 3))` -- tick must reach 3 again and
// again -- and it therefore FAILS.  Obliging the property only in the
// first cycle would instead ask for "tick == 3 at some point", which
// is true, and report no violation.
module eventually_not_recurring (
    input logic clk,
    input logic rst
);

  logic [3:0] tick;

  always_ff @(posedge clk) begin
    if (rst) tick <= 4'd0;
    else if (tick != 4'hF) tick <= tick + 4'd1;
  end

  assert property (@(posedge clk) s_eventually (tick == 4'd3));

endmodule
