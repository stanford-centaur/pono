// Out-of-scope: functional coverage (`covergroup`/`coverpoint`).
module covergroup_decl (input logic clk, input logic rst, input logic [3:0] a);
  covergroup cg @(posedge clk);
    coverpoint a;
  endgroup
  cg cg_inst = new();
  assert property (@(posedge clk) 1'b1);
endmodule
