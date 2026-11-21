// The inline twin of named_property_expr_arg.sv.
module named_property_expr_arg_inline (
    input logic clk,
    input logic req,
    input logic ack
);

  assert property (@(posedge clk) (req && !ack) |-> ##1 ack);

endmodule
