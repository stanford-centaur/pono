// An explicit port bound to a bit range rather than a whole signal:
// `output .o(w[3:0])` carries four of `w`'s eight bits. An output
// alias already maps a port range onto a range of a target, so the
// shape has a representation; what it needs is for the binding to be
// read from the expression.
//
// The child drives its port, so the low half of `w` follows `d` and
// the high half stays as the design left it.
module explicit_port_partial_child (
    input  logic [3:0] d,
    output logic [3:0] o
);
  assign o = d;
endmodule

module explicit_port_partial (
    input logic clk,
    input logic [7:0] d,
    output .o(w[3:0])
);

  logic [7:0] w;
  assign w[7:4] = 4'hA;

  explicit_port_partial_child u (.d(d[3:0]), .o(w[3:0]));

  assert property (@(posedge clk) w[3:0] == d[3:0] && w[7:4] == 4'hA);

endmodule
