// A multi-element sequence used directly as a property
// (`assert property (a ##1 b);`), as opposed to as the antecedent of
// `|->`/`|=>` (which assertion_expr_to_bool() does handle via the
// bounded sequence matcher -- see e.g. SeqIntersect and friends).
// Per the LRM this means "the sequence must eventually match", an
// inherent liveness obligation ltl_to_sat()'s SequenceConcat case has
// no gadget for.
module bare_sequence_concat_property (
    input logic clk,
    input logic a,
    input logic b
);
  assert property (@(posedge clk) a ##1 b);
endmodule
