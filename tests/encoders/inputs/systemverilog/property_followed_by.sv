// Overlapped followed-by (`#-#`, `AssertionExprKind::Binary` with
// `BinaryAssertionOperator::OverlappedFollowedBy`) as a top-level
// property connective -- same unmodeled-multi-cycle-operator family
// as intersect/within/throughout, see property_intersect_toplevel.sv.
module property_followed_by (
    input logic clk,
    input logic a,
    input logic b
);
  assert property (@(posedge clk) a #-# b);
endmodule
