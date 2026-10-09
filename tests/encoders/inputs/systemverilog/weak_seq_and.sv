// weak() of `and` over operands of different lengths. The operands
// start together and the composite ends where the later one ends
// (LRM 16.9.5), so `(x ##1 y) and (x ##2 z)` is the same match set
// as the chain `x ##1 y ##1 z` in weak_seq_and_chain.sv -- the two
// must be refuted at the same depth, which is the span of the
// *longer* operand rather than either the shorter one's or the sum.
module weak_seq_and (
    input logic clk,
    input logic x,
    input logic y,
    input logic z
);

  assert property (@(posedge clk) weak((x ##1 y) and (x ##2 z)));

endmodule
