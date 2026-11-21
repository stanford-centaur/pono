// An explicit (named) output port: `.o(r)` gives the port the
// outward name `o` while binding it to the internal signal `r`.
// slang leaves PortSymbol::internalSymbol null for this form and
// puts the binding in an expression instead, so reading that field
// directly found nothing and skipped the connection in silence --
// leaving the parent's `q` with no driver, which showed up only if
// something happened to read it.
//
// Paired with explicit_port_out_plain.sv, the same design written
// with an ordinary port. The two must prove identically.
module explicit_port_out_sub (
    input logic clk,
    input logic d,
    output .o(r)
);
  logic r;
  always_ff @(posedge clk) r <= d;
endmodule

module explicit_port_out (
    input logic clk,
    input logic d
);

  logic q;
  logic started;
  initial started = 1'b0;
  always_ff @(posedge clk) started <= 1'b1;

  explicit_port_out_sub u (
      .clk(clk),
      .d  (d),
      .o  (q)
  );

  assert property (@(posedge clk) started |-> q == $past(d));

endmodule
