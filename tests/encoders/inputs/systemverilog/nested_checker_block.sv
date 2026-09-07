// A checker instantiated inside a non-top module, rather than
// directly in the top module the way checker_block.sv does --
// exercises process_instance()'s own CheckerInstance dispatch.
checker my_checker (input logic clk, input logic a);
  assert property (@(posedge clk) a);
endchecker

module child (input logic clk, input logic a);
  my_checker chk (clk, a);
endmodule

module nested_checker_block (input logic clk, input logic a);
  child c (clk, a);
endmodule
