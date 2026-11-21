// A named sequence referenced from inside another named assertion.
// The two property-level callers resolved such a reference already,
// but reaching one as a *sequence* operand -- here as an
// implication's antecedent -- landed on expr_to_bool(), which knows
// nothing of assertion instances and reported "unsupported
// expression kind 39": a raw slang enum index with no source
// location, which is what made it hard to place.
//
// Paired with named_assertion_nested_inlined.sv, the same property
// written out. Same verdict at the same depth, or the reference is
// resolving to something other than its body.
module named_assertion_nested (
    input logic clk,
    input logic x,
    input logic y
);

  sequence inner(s);
    s ##1 s;
  endsequence

  property outer(p, q);
    inner(p) |=> q;
  endproperty

  assert property (@(posedge clk) outer(x, y));

endmodule
