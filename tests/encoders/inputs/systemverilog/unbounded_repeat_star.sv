// `a[*]` (`[*0:$]`) as the whole antecedent of `|->`. It matches
// emptily, and an empty match occupies no cycles, so it matches at
// *every* cycle: the implication is unconditional and the `a` plays
// no part. Degenerate, but well defined -- the property is just its
// consequent -- so this encodes, with a warning, rather than being
// refused. (Which cycle an empty match is said to end on does not
// change the answer: an assertion is checked at every cycle, so
// "c at k" and "c at k-1" for all k are the same constraint.)
//
// `c` is false exactly at cnt == 5, so this is refuted there --
// the same cycle as a bare `assert property (c)`, which is what
// unbounded_repeat_star_equiv.sv pins.
module unbounded_repeat_star (
    input logic clk,
    input logic rst
);

  logic [3:0] cnt;
  logic a, c;

  always_ff @(posedge clk) cnt <= rst ? 4'd0 : cnt + 4'd1;
  assign a = !rst && (cnt == 4'd3);
  assign c = !rst && (cnt != 4'd5);

  assert property (@(posedge clk) (!rst) |-> (a[*] |-> c));

endmodule
