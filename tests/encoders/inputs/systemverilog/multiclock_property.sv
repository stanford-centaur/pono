// A multiclock property (`@(posedge clk1) a ##1 @(posedge clk2) b`).
// The Clocking-wrapper handling in the ConcurrentAssertion/
// assertion_expr_to_bool code strips *every* clocking wrapper it
// finds unconditionally, without checking whether the newly-named
// clock differs from the outer one -- so a genuine clock change may
// be silently collapsed into a single-clock (mis-)encoding rather
// than either being honored or cleanly rejected.  Verified
// empirically which of those three outcomes actually happens.
module multiclock_property (
    input logic clk1,
    input logic clk2,
    input logic rst,
    input logic a,
    input logic b
);

  assert property (@(posedge clk1) a ##1 @(posedge clk2) b |-> 1'b0);

endmodule
