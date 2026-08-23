// `cover sequence (...)`, distinct from `cover property (...)`
// (already supported via reachability duality) -- both share
// AssertionKind::CoverProperty vs. AssertionKind::CoverSequence, and
// the ConcurrentAssertion statement handler only ever checked for the
// latter, silently dropping this construct entirely (no log, no
// throw, no encoding effect).
module cover_sequence (input logic clk, input logic a, input logic b);

  cover sequence (@(posedge clk) a ##1 b);

endmodule
