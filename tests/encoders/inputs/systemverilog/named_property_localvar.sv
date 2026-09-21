// A local variable inside a named property needs a binding
// environment of its own, which expanding the body does not supply,
// so this still has to be rejected.
module named_property_localvar (
    input logic clk,
    input logic [7:0] d,
    input logic v
);

  property saved;
    logic [7:0] t;
    (v, t = d) |-> ##1 (d == t);
  endproperty

  assert property (@(posedge clk) saved);

endmodule
