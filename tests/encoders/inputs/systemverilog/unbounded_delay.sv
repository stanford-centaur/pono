// `req |-> ##[1:$] gnt`: the canonical liveness property. An
// unbounded wait has no finite window for the bounded sequence
// matcher to span, so the property goes to the LTL tableau, where
// the delay fixes where the wait starts and F carries it the rest of
// the way.
//
// The grant always arrives the cycle after the request, so every
// eventuality is discharged and no fair lasso violates this.
// unbounded_delay_fails.sv is the paired refutation: an F that never
// discharges its eventuality would pass here too.
module unbounded_delay (
    input logic clk,
    input logic req
);

  logic gnt;

  always_ff @(posedge clk) gnt <= req;

  assert property (@(posedge clk) req |-> ##[1:$] gnt);

endmodule
