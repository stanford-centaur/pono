// How far an abort reaches. `a |=> b` is checked one cycle after
// `a`, so the attempt beginning at cnt==0 fails at cnt==1; `abrt`
// pulses only at cnt==0, the attempt's *first* cycle. The abort
// covers the whole attempt, so that failure is waived.
//
// This is the case that catches a window derived at the wrong
// level. An implication gates itself across its own
// antecedent-to-consequent span and then reports a span of zero,
// having nothing left for its caller to do -- so an abort that
// built its own window from that report would waive only the cycle
// the check lands on, and this design would be refuted.
module property_accept_on_window (input logic clk);

  logic [3:0] cnt;
  initial cnt = 4'd0;
  always @(posedge clk) cnt <= cnt + 4'd1;

  logic abrt, a, b;
  assign abrt = (cnt == 4'd0);
  assign a   = (cnt == 4'd0);
  assign b   = 1'b0;

  assert property (@(posedge clk) accept_on (abrt) (a |=> b));

endmodule
