// The same sequence and the same design as
// bare_sequence_concat_property.sv, with the `weak` qualifier
// `assert property` supplies for it written out. LRM 16.13.1 makes
// the two spellings the same property, so they must be refuted at
// the same depth; the encoder used to prove this one and refute the
// other.
//
// `x` is tied low, so no attempt can ever match from any tick, and
// no assignment of the free `y` can rescue it.
module weak_seq_never_matches (
    input logic clk,
    input logic y
);

  logic x;
  assign x = 1'b0;

  assert property (@(posedge clk) weak(x ##1 y));

endmodule
