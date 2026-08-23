// Sequence `throughout`: per the LRM, `expr throughout seq` requires
// the (single-cycle, boolean) expr to hold at *every* cycle spanned by
// seq's match. For `b ##1 c` (spanning [i, i+1]), that reduces to:
// a@i AND a@i+1 AND b@i AND c@i+1. Wrapped in `|-> 1'b0` (the same
// idiom first_match_seq.sv/seq_intersect.sv/seq_within.sv use) so the
// property is violated as soon as the composite sequence matches at
// all. `a`/`b`/`c` free: earliest match is a@1, b@1, a@2, c@2, so the
// property is violated at cycle 2. Handled by offsets_ending_now()'s
// Throughout case in sva.cpp, which ANDs the plain boolean over seq's
// own completion window.
module seq_throughout (input logic clk, input logic rst, input logic a, input logic b, input logic c);

  assert property (@(posedge clk) (a throughout (b ##1 c)) |-> 1'b0);

endmodule
