// Sequence `throughout`: the same silent-drop AssertionExprKind gap
// as seq_intersect.sv, for a different sequence composition operator.
module seq_throughout (input logic clk, input logic rst, input logic a, input logic b, input logic c);

  assert property (@(posedge clk) a throughout (b ##1 c));

endmodule
