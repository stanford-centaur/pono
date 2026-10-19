// `join_none`: the parent continues at once and the spawned
// processes do not start until it blocks (LRM 9.3.2, Table 9-1). So
// the writes below have not happened when the parent moves on --
// running them there would be executing code the language says has
// not begun.
module fork_join_none (
    input logic clk
);

  logic [3:0] x, y;

  initial begin
    x = 4'd0;
    y = 4'd0;
    fork
      x = 4'd5;
      y = 4'd9;
    join_none
  end

  assert property (@(posedge clk) x == 4'd5 && y == 4'd9);

endmodule
