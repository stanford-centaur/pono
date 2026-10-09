// SVA `and` operator: both operands must hold this cycle.
// Equivalent to logical conjunction for simple boolean operands.
// Free input `data` is the antecedent of both halves; BMC picks
// data so that at least one half fails to falsify.
module binary_and (input logic clk, input logic [3:0] data);
  assert property (@(posedge clk) (data[0] == 1'b1) and (data[1] == 1'b0));
endmodule
