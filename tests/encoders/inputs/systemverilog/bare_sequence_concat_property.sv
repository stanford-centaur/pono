// A multi-element sequence used directly as a property. An
// evaluation attempt is a search for a match *beginning at a
// particular clock tick* (LRM 16.10), so this is not an eventuality:
// at each tick the attempt starting there must match, and `a` being
// low at that tick already fails it.
//
// Weak and strong differ only over whether the remaining ticks
// exist. 16.13.1's worked example is explicit that weak still
// requires the sequence's first element at the attempt's own tick --
// and `assert property` takes the weak reading by default, which
// changed from strong in IEEE 1800-2009. On the infinite traces Pono
// reasons about the later ticks always exist, so the two coincide
// and this is one per-cycle safety check.
//
// With `a` free it is refuted as soon as a tick has `a` low.
module bare_sequence_concat_property (
    input logic clk,
    input logic a,
    input logic b
);
  assert property (@(posedge clk) a ##1 b);
endmodule
