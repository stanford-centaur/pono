// The discriminating negative for explicit_port_in.sv: the child
// inverts what it was given, so the same assertion must now be
// refuted. Without this, an explicit input port left unconstrained
// would still pass the positive test -- a free variable satisfies
// `q == d` on some trace.
module explicit_port_in_fails_sub (
    input logic clk,
    input .i(a),
    output logic o
);
  logic a;
  assign o = ~a;
endmodule

module explicit_port_in_fails (
    input logic clk,
    input logic d
);

  logic q;

  explicit_port_in_fails_sub u (
      .clk(clk),
      .i  (d),
      .o  (q)
  );

  assert property (@(posedge clk) q == d);

endmodule
