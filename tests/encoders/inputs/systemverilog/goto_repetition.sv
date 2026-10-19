// `gnt[->1]`: the next occurrence of gnt. Counting occurrences that
// need not be adjacent spans no finite window, so the bounded
// matcher declines and the tableau reads it as the eventuality it
// is. The grant always arrives the cycle after the request, so every
// one of them is reached.
//
// `[=1]` says the same thing here: its match may run on past the
// occurrence, but the earliest one ends there, and a consequent only
// has to match somewhere.
module goto_repetition (
    input logic clk,
    input logic req
);

  logic gnt;

  always_ff @(posedge clk) gnt <= req;

  assert property (@(posedge clk) req |-> gnt[->1] and gnt[=1]);

endmodule
