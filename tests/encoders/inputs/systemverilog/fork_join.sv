// Out-of-scope: `fork`/`join` concurrent processes.
module fork_join (input logic clk, input logic rst);
  initial begin
    fork
      #1;
      #2;
    join
  end
  assert property (@(posedge clk) 1'b1);
endmodule
