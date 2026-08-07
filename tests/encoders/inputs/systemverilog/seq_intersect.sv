// Sequence `intersect`: both operand sequences must match with the
// *same* length.  Absent from the AssertionExprKind switches in
// assertion_expr_to_bool()/ltl_to_sat() (confirmed by inspection --
// only Clocking/StrongWeak/Simple/SequenceConcat/Unary/Binary are
// handled), so the default case returns a null Term and the caller
// silently drops the whole assertion rather than throwing.
module seq_intersect (input logic clk, input logic rst, input logic a, input logic b);

  assert property (@(posedge clk) (a ##1 b) intersect (a ##[1:2] b));

endmodule
