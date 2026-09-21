// named_assertion_nested.sv with both references written out.
module named_assertion_nested_inlined (
    input logic clk,
    input logic x,
    input logic y
);

  assert property (@(posedge clk) (x ##1 x) |=> y);

endmodule
