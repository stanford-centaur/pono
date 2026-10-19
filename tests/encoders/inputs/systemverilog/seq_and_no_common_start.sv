// The negative of seq_and_common_start.sv, and the one that catches
// a model too permissive rather than too strict: `a ##1 b` can only
// match starting at cnt==2 and `c ##3 d` only starting at cnt==0, so
// there is no cycle both operands start at and `and` never matches.
// Proved rather than merely left unrefuted -- a bounded run cannot
// tell "no match" from "no match yet".
module seq_and_no_common_start (input logic clk);

  logic [3:0] cnt;
  initial cnt = 4'd0;
  always @(posedge clk) cnt <= cnt + 4'd1;

  logic a, b, c, d;
  assign a = (cnt == 4'd2);
  assign b = (cnt == 4'd3);
  assign c = (cnt == 4'd0);
  assign d = (cnt == 4'd3);

  assert property (@(posedge clk) ((a ##1 b) and (c ##3 d)) |-> 1'b0);

endmodule
