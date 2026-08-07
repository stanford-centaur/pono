// The postfix increment statement (`cnt++;`) as the sole next-state
// update for a register -- distinct from the `cnt <= cnt + 1`
// spelling used elsewhere in this suite.  GAP (confirmed empirically,
// not just assumed): `cnt++;` is a bare ExpressionStatement whose
// expression kind is UnaryOp/Postincrement, not Assignment, and
// process_statement's ExpressionStatement case only recognizes
// ExpressionKind::Assignment -- so this line is silently a no-op and
// cnt gets stuck at its reset value (0) forever, rather than throwing.
// The intended semantics (equivalent to `cnt <= cnt + 1`) would reach
// cnt == 3 at cycle 4; instead the property holds vacuously at every
// bound, which is exactly the dangerous "elaborates fine, drops
// semantics" failure mode called out in the plan's Context section.
module inc_dec (input logic clk, input logic rst);

  logic [3:0] cnt;

  always_ff @(posedge clk) begin
    if (rst) cnt <= 4'd0;
    else cnt++;
  end

  assert property (@(posedge clk) cnt != 4'd3);

endmodule
