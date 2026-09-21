// A named *sequence* with arguments, used directly as a property.
// A bare bounded sequence is a per-cycle safety check, and this
// takes that path exactly as its inline twin does -- it used to
// reach the LTL tableau instead, because the reference was opaque to
// the sequence matcher. Paired with named_sequence_args_inline.sv,
// which it must match down to the refutation depth.
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
