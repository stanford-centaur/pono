// `a ##[1:$] b` as an antecedent. The gap between the two elements
// has no bound, so no window of offsets spans the match -- but where
// the match *ends* is still a definite cycle, and that is all an
// antecedent needs. One latch remembers that `a` happened; `b`
// supplies the end.
//
// So the antecedent fires exactly when `b` holds and `a` held at
// some earlier cycle, which `seen_a` states by hand. The stimulus is
// gated on reset because the encoder's latch starts at time 0 while
// this fixture's is cleared by `rst` -- ungated, the two would
// disagree about a free input firing during reset.
//
// Paired with unbounded_delay_antecedent_fails.sv.
module unbounded_delay_antecedent (
    input logic clk,
    input logic rst,
    input logic ai,
    input logic bi
);

  logic a, b, seen_a;

  assign a = ai && !rst;
  assign b = bi && !rst;
  always_ff @(posedge clk) seen_a <= rst ? 1'b0 : (seen_a || a);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> ((a ##[1:$] b) |-> (b && seen_a)));

endmodule
