// A repetition applied to the reference itself. `inner(x)[*2]`
// repeats a whole sequence, which is a convolution of its offsets
// with themselves -- not the per-cycle repetition a Boolean operand
// gets. Refused rather than quietly given the wrong one, which is
// the reading a bare `resolve and recurse` would have fallen into.
module named_assertion_reference_repetition (
    input logic clk,
    input logic x,
    input logic y
);

  sequence inner(s);
    s ##1 s;
  endsequence

  assert property (@(posedge clk) inner(x)[*2] |=> y);

endmodule
