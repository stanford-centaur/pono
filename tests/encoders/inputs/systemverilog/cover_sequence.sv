// `cover sequence (...)` is treated the same as `cover property (...)`
// (reachability duality). `a ##1 b` is multi-cycle but not an
// eventuality: an attempt searches for a match beginning at its own
// tick, so the goal is a per-cycle check like any other.
module cover_sequence (input logic clk, input logic a, input logic b);

  cover sequence (@(posedge clk) a ##1 b);

endmodule
