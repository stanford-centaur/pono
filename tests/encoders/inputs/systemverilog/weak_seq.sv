// `weak(seq)` where the continuation after the leading element is an
// unconditional truth, so the only thing that can fail the attempt is
// the leading element itself.
//
// An attempt begins at *every* tick (LRM 16.10), not only where the
// leading element happens to hold, and `weak(R)` fails as soon as a
// finite prefix witnesses inability to match (16.13.1, F.5.3.1).
// With `a` low at a tick, that tick's own attempt is already dead --
// so this is exactly as strong as `a` itself, and `a` free means it
// is refuted as soon as a tick has `a` low. Same reading, and the
// same depth, as the unwrapped form in
// bare_sequence_concat_property.sv, which `assert property` makes
// weak anyway.
module weak_seq (input logic clk, input logic rst, input logic a);

  assert property (@(posedge clk) weak(a ##1 1'b1));

endmodule
