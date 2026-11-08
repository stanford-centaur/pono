// A fork whose branches carry timing of their own, which is what
// parallel blocks are for: "the timing controls in a fork-join block
// do not have to be ordered sequentially in time" (LRM 9.3.2). There
// is no time between clock edges here for the delays to land in, and
// no scheduler to run the branches against each other.
module fork_join_timed (
    input logic clk
);

  logic [3:0] x, y;

  initial begin
    x = 4'd0;
    y = 4'd0;
    fork
      #2 x = 4'd5;
      #1 y = 4'd9;
    join
  end

  assert property (@(posedge clk) x == 4'd5 && y == 4'd9);

endmodule
