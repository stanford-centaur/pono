// SVA `iff` operator on plain boolean operands: equivalent to
// logical equivalence. Mirrors binary_and.sv's structure.
module binary_iff (input logic clk, input logic [3:0] data);
  assert property (@(posedge clk) (data[0] == 1'b1) iff (data[1] == 1'b0));
endmodule
