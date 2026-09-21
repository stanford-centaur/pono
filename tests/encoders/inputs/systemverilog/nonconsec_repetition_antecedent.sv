// `b[=n]` as an antecedent, which differs from `b[->n]` exactly in
// that it need not end on an occurrence: the same counter, without
// the requirement that `b` holds now. So this must be refuted --
// `b` is low at cnt == 6, one cycle after the 2nd occurrence, and
// the antecedent still fires there.
module nonconsec_repetition_antecedent (
    input logic clk,
    input logic rst
);

  logic [3:0] cnt;
  logic b;

  always_ff @(posedge clk) cnt <= rst ? 4'd0 : cnt + 4'd1;
  assign b = !rst && ((cnt == 4'd2) || (cnt == 4'd5) || (cnt == 4'd9));

  assert property (@(posedge clk) (!rst) |-> (b[=2] |-> b));

endmodule
