// Two levels of nesting, with the inner reference sitting inside a
// sequence concatenation rather than alone -- so the resolution has
// to happen where offsets_ending_now() recurses into a concat
// element, not only at the top of a sequence.
module named_assertion_nested_deep (
    input logic clk,
    input logic x,
    input logic y
);

  sequence lvl1(s);
    s ##1 s;
  endsequence

  sequence lvl2(s);
    lvl1(s) ##1 s;
  endsequence

  property top(p, q);
    lvl2(p) |=> q;
  endproperty

  assert property (@(posedge clk) top(x, y));

endmodule
