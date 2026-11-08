// A named property taking arguments -- the parameterized-property
// idiom most real SVA is written in. Slang expands the body with the
// actuals already substituted, so the encoder recurses into it the
// same way it does for an argument-free reference.
//
// Paired with named_property_args_inline.sv, which spells the same
// property out: both must be refuted at the same depth, which is
// what shows the arguments arrived where they belong rather than the
// reference merely encoding.
module named_property_args (
    input logic clk,
    input logic req,
    input logic ack
);

  property handshake(a, b);
    a |-> ##1 b;
  endproperty

  assert property (@(posedge clk) handshake(req, ack));

endmodule
