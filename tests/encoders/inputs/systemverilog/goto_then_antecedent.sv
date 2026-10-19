// A counted repetition followed by another element, as an
// antecedent. The count yields a Boolean "a match ends now", and
// `##1 b` after it needs only that Boolean delayed a cycle -- so
// the whole thing fires when `b` holds and `a[->2]` matched on the
// previous cycle. `goto2` states that by hand.
//
// A leading `##[1:2]` is also covered here: a delay before the
// count moves where the attempt started, and a match whose start is
// already an unbounded distance back is indifferent to that.
module goto_then_antecedent (
    input logic clk,
    input logic rst,
    input logic ai,
    input logic bi
);

  logic a, b, goto2;
  logic [1:0] na;

  assign a = ai && !rst;
  assign b = bi && !rst;
  always_ff @(posedge clk)
    na <= rst ? 2'd0 : ((a && na < 2'd2) ? na + 2'd1 : na);
  assign goto2 = a && na >= 2'd1;

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> (((a[->2] ##1 b) |-> (b && $past(goto2)))
                        and ((##[1:2] a[->2]) |-> goto2)));

endmodule
