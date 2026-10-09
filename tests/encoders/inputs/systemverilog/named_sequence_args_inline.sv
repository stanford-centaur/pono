// The inline twin of named_sequence_args.sv.
module named_sequence_args_inline (
    input logic clk,
    input logic req,
    input logic ack
);

  assert property (@(posedge clk) req ##1 ack);

endmodule
