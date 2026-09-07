// Out-of-scope: a procedural `wait` statement.
module wait_stmt (input logic clk, input logic rst, input logic go);
  initial begin
    wait (go);
  end
  assert property (@(posedge clk) 1'b1);
endmodule
