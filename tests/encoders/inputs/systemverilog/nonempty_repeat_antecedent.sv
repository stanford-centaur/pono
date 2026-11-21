// The control for unbounded_repeat_star.sv: the same shape with a
// *non-zero* lower bound, so there is no empty match and the
// antecedent still depends on `a`. `a` is true only at cnt == 3 and
// `c` is false only at cnt == 5, so no match ever meets a false
// consequent and this is not refuted -- which it would be if the
// collapse applied to every repetition rather than the empty one.
module nonempty_repeat_antecedent (
    input logic clk,
    input logic rst
);

  logic [3:0] cnt;
  logic a, c;

  always_ff @(posedge clk) cnt <= rst ? 4'd0 : cnt + 4'd1;
  assign a = !rst && (cnt == 4'd3);
  assign c = !rst && (cnt != 4'd5);

  assert property (@(posedge clk) (!rst) |-> (a[*1:2] |-> c));

endmodule
