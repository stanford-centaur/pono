// An in-property `if`/`else` (`AssertionExprKind::Conditional`) as a
// top-level property connective -- distinct from a plain procedural
// `if` inside an always block, this chooses which sub-property to
// check based on `sel`, which no gadget in this tableau builds.
module property_conditional (
    input logic clk,
    input logic sel,
    input logic a,
    input logic b
);
  assert property (@(posedge clk) if (sel) a else b);
endmodule
