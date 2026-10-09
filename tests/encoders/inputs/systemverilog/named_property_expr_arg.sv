// An argument that is an expression rather than a bare signal, so
// the substitution has to carry a whole subexpression into the body.
// Paired with named_property_expr_arg_inline.sv.
module named_property_expr_arg (
    input logic clk,
    input logic req,
    input logic ack
);

  property handshake(a, b);
    a |-> ##1 b;
  endproperty

  assert property (@(posedge clk) handshake(req && !ack, ack));

endmodule
