// `disable iff` on a property that goes through the LTL tableau
// rather than the safety fast path (`nexttime` is temporal, so
// assertion_expr_to_bool() declines it).
//
// `started` is cleared while rst is high, so it is still low one cycle
// after the reset cycle.  `nexttime started` therefore fails when
// evaluated at the reset cycle, and holds at every later cycle.  The
// `disable iff (rst)` exempts exactly that cycle, so the assertion
// holds.  Paired with disable_iff_temporal_fails.sv, which is the same
// design and property with the exemption removed.
module disable_iff_temporal (
    input logic clk,
    input logic rst
);

  logic started;

  always_ff @(posedge clk) begin
    if (rst) started <= 1'b0;
    else started <= 1'b1;
  end

  assert property (@(posedge clk) disable iff (rst) nexttime started);

endmodule
