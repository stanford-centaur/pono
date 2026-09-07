// Sequence `intersect`: both operand sequences must match with the
// *same* length.  `a ##1 b` only ever matches with length 1 (a@i,
// b@i+1); `a ##[1:2] b` matches with length 1 (a@i, b@i+1) or length 2
// (a@i, b@i+2) -- the only length both can share is 1, so the intersect
// is logically equivalent to plain `a ##1 b`: it matches iff a holds
// and b holds exactly one cycle later.  Wrapped in `|-> 1'b0` (the same
// idiom first_match_seq.sv uses) so the property is violated as soon as
// the sequence matches at all, sidestepping the separate (and more
// LRM-ambiguous) question of what a *bare* sequence used directly as a
// property means. `a`/`b` free: earliest match is a@1, b@2, so the
// property is violated at cycle 2. Handled by offsets_ending_now()'s
// Intersect case in sva.cpp, which ANDs the two operands' offset
// vectors entry-by-entry.
module seq_intersect (input logic clk, input logic rst, input logic a, input logic b);

  assert property (@(posedge clk) ((a ##1 b) intersect (a ##[1:2] b)) |-> 1'b0);

endmodule
