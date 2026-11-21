// `wait fork;` blocks the calling process until every one of its
// child subprocesses has finished (LRM 9.6.1). It is a statement
// about processes, and there are no processes here to wait on.
module wait_fork (
    input logic clk
);

  logic [3:0] x;

  initial begin
    x = 4'd0;
    wait fork;
    x = 4'd5;
  end

  assert property (@(posedge clk) x == 4'd5);

endmodule
