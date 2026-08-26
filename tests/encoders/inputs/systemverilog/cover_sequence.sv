// `cover sequence (...)` is treated the same as `cover property (...)`
// (reachability duality). Since `a ##1 b` is a genuinely multi-cycle
// sequence, it hits the temporal/sequence-shaped cover-goal throw.
module cover_sequence (input logic clk, input logic a, input logic b);

  cover sequence (@(posedge clk) a ##1 b);

endmodule
