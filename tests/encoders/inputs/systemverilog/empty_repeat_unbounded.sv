// `a[*0:$]` needs nothing beyond what `a[*1:$]` already does once the
// empty match is reported separately: a run of at least one `a`
// ending now ends with exactly one, so every non-empty alternative
// contributes the same offset. Every alternative here still ends on
// `b`, so the whole sequence completes exactly when `b` holds.
//
// `a[*]` on its own remains rejected -- see unbounded_repeat_star.sv.
module empty_repeat_unbounded (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b
);

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> ((a[*0:$] ##1 b) |-> b));

endmodule
