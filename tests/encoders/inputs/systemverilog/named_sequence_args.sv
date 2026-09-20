// A named *sequence* with arguments, used directly as a property.
// That is a bare sequence, so it goes to the LTL tableau rather than
// the safety path -- the argument substitution has to hold up there
// too. Paired with named_sequence_args_inline.sv.
module named_sequence_args (
    input logic clk,
    input logic req,
    input logic ack
);

  sequence handshake(a, b);
    a ##1 b;
  endsequence

  assert property (@(posedge clk) handshake(req, ack));

endmodule
