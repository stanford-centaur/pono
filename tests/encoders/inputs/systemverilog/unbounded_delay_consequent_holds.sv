// unbounded_delay_consequent.sv with `b` available, so the
// eventuality can always be discharged and no fair lasso violates
// it. Without this pair the refutation there could just as well be
// the encoder making the consequent unsatisfiable outright.
module unbounded_delay_consequent_holds (
    input logic clk,
    input logic rst
);

  logic a, b, c;

  assign a = !rst;
  assign b = !rst;
  assign c = !rst;

  assert property (@(posedge clk) c |-> (a ##[1:$] b));

endmodule
