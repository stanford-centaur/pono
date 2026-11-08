// The same two operands as seq_and_no_common_start.sv under `or`
// instead. A union of match sets asks nothing of where the other
// operand starts, so the pair that `and` rejects still matches here
// -- at cnt==3, where both alternatives happen to end.
module seq_or_independent_matches (input logic clk);

  logic [3:0] cnt;
  initial cnt = 4'd0;
  always @(posedge clk) cnt <= cnt + 4'd1;

  logic a, b, c, d;
  assign a = (cnt == 4'd2);
  assign b = (cnt == 4'd3);
  assign c = (cnt == 4'd0);
  assign d = (cnt == 4'd3);

  assert property (@(posedge clk) ((a ##1 b) or (c ##3 d)) |-> 1'b0);

endmodule
