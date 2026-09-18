// `accept_on` (`AssertionExprKind::Abort`, action=Accept) as a
// top-level property connective. accept_on/reject_on/sync_accept_on/
// sync_reject_on all share this AssertionExprKind, distinguished only
// by the Action enum/isSync flag -- none are modeled.
module property_accept_on (
    input logic clk,
    input logic rst,
    input logic a
);
  assert property (@(posedge clk) accept_on (rst) a);
endmodule
