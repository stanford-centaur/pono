// A repetition applied to a named sequence reference. `inner(x)[*2]`
// repeats a whole sequence, which is its offsets convolved with
// themselves -- the same composition `(seq)[*n]` written inline gets,
// rather than the per-cycle repetition a Boolean operand takes.
//
// `inner(x)` is `x ##1 x`, so `inner(x)[*2]` is four consecutive
// cycles of `x`, and the consequent is checked at the fifth.
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
