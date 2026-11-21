// `fork ... join` spawns each statement as a process and waits for
// all of them (LRM 9.3.2). Concurrent processes with timing of their
// own are not what a cycle-based transition system describes, and
// running the branches one after another instead would settle a race
// the language leaves open -- two branches writing one variable have
// no defined winner.
//
// The branches here write different variables, so this particular
// fork happens to mean the same as `begin ... end`. That spelling is
// available and says so.
module fork_join (
    input logic clk
);

  logic [3:0] x, y;

  initial begin
    x = 4'd0;
    y = 4'd0;
    fork
      x = 4'd5;
      y = 4'd9;
    join
  end

  assert property (@(posedge clk) x == 4'd5 && y == 4'd9);

endmodule
