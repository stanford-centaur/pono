// The inline twin of named_property_args.sv.
module named_property_args_inline (
    input logic clk,
    input logic req,
    input logic ack
);

  assert property (@(posedge clk) req |-> ##1 ack);

endmodule
