// A minimum above one, so the wait starts two cycles out rather
// than at the next. The grant arrives exactly one cycle after the
// request and never again, which is too early to satisfy
// `##[2:$]` -- so this must be refuted, and would not be if the
// minimum were dropped on the way to the tableau.
module unbounded_delay_min (
    input logic clk,
    input logic req
);

  logic gnt;

  always_ff @(posedge clk) gnt <= req && !gnt;

  assert property (@(posedge clk) req |-> ##[2:$] gnt);

endmodule
