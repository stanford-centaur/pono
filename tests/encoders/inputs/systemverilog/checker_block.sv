// Out-of-scope: a `checker` block instantiated like a module.
checker my_checker (input logic clk, input logic a);
  assert property (@(posedge clk) a);
endchecker

module checker_block (input logic clk, input logic rst, input logic a);
  my_checker chk (clk, a);
endmodule
