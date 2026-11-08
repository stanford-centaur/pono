// Each signal is true at exactly one cycle, so every operand has
// exactly one possible match and where it starts is pinned. Both
// `a ##1 b` and `c ##3 d` start at cnt==0 here, so `and` matches,
// ending with the later of the two -- at cnt==3. A consequent of
// 1'b0 turns "the antecedent matched" into the counterexample, so
// the refutation depth *is* the composite's end cycle. A model that
// merged the operands without bringing the shorter one back to the
// common start would look for `a ##1 b` ending at cnt==3 instead and
// find nothing.
module seq_and_common_start (input logic clk);

  logic [3:0] cnt;
  initial cnt = 4'd0;
  always @(posedge clk) cnt <= cnt + 4'd1;

  logic a, b, c, d;
  assign a = (cnt == 4'd0);
  assign b = (cnt == 4'd1);
  assign c = (cnt == 4'd0);
  assign d = (cnt == 4'd3);

  assert property (@(posedge clk) ((a ##1 b) and (c ##3 d)) |-> 1'b0);

endmodule
