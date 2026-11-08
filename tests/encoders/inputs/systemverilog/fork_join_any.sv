// `join_any`: the parent resumes as soon as *one* branch finishes,
// with the others still running (LRM 9.3.2, Table 9-1). Carrying a
// process that has started and not finished past the point where the
// parent continues needs a notion of a process, which a transition
// system does not have.
module fork_join_any (
    input logic clk
);

  logic [3:0] x, y;

  initial begin
    x = 4'd0;
    y = 4'd0;
    fork
      x = 4'd5;
      y = 4'd9;
    join_any
  end

  assert property (@(posedge clk) x == 4'd5 && y == 4'd9);

endmodule
