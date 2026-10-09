// An unbounded inter-element delay in a *consequent*. There the
// sequence is an eventuality, not a safety obligation, so the
// bounded-span matcher declines and the tableau takes over: the
// obligation is that a match ends at some point, and "a match ends
// at this cycle" is what the encoder supplies however far back the
// match began.
//
// `b` never holds, so no match can ever complete and the
// eventuality is unsatisfiable -- a fair lasso exists and this is
// refuted. unbounded_delay_consequent_holds.sv is the same property
// where `b` is available, which is not.
module unbounded_delay_consequent (
    input logic clk,
    input logic rst
);

  logic a, b, c;

  assign a = !rst;
  assign b = 1'b0;
  assign c = !rst;

  assert property (@(posedge clk) c |-> (a ##[1:$] b));

endmodule
