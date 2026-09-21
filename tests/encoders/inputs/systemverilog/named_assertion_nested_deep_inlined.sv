// named_assertion_nested_deep.sv with every reference written out.
module named_assertion_nested_deep_inlined (
    input logic clk,
    input logic x,
    input logic y
);

  assert property (@(posedge clk) (x ##1 x ##1 x) |=> y);

endmodule
