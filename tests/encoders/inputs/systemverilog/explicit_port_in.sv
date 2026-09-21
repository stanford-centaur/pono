// An explicit *input* port. The failure here was quieter than the
// output one and pointed the other way: the child's `a` was left an
// unconstrained variable rather than the parent's `d`, so a design
// that trivially holds was refuted.
module explicit_port_in_sub (
    input logic clk,
    input .i(a),
    output logic o
);
  logic a;
  assign o = a;
endmodule

module explicit_port_in (
    input logic clk,
    input logic d
);

  logic q;

  explicit_port_in_sub u (
      .clk(clk),
      .i  (d),
      .o  (q)
  );

  assert property (@(posedge clk) q == d);

endmodule
