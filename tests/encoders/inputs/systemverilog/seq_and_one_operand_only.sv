// `and` requires *both* operands. `c` is true only at cnt==0 and
// `d` only at cnt==9, which is nine cycles later rather than the
// three `##3` allows, so `c ##3 d` has no match at all and the
// composite has none either -- however well `a ##1 b` does on its
// own, which here is once every time the counter wraps.
module seq_and_one_operand_only (input logic clk);

  logic [3:0] cnt;
  initial cnt = 4'd0;
  always @(posedge clk) cnt <= cnt + 4'd1;

  logic a, b, c, d;
  assign a = (cnt == 4'd0);
  assign b = (cnt == 4'd1);
  assign c = (cnt == 4'd0);
  assign d = (cnt == 4'd9);

  assert property (@(posedge clk) ((a ##1 b) and (c ##3 d)) |-> 1'b0);

endmodule
