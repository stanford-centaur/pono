// A checker that declares its own local variable, rather than only
// referencing its formal ports (which slang substitutes through to
// the caller's own symbols, needing no declaration on this side at
// all -- see checker_block.sv). `local_var` is genuine new state this
// encoder would need its own pre-scan/declare pass to handle.
checker my_checker (input logic clk, input logic a);
  logic local_var;
  assert property (@(posedge clk) a);
endchecker

module checker_local_variable (input logic clk, input logic a);
  my_checker chk (clk, a);
endmodule
