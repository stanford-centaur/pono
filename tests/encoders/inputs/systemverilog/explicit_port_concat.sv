// An explicit port bound to a concatenation of several signals,
// `output .o({x, y})`. An ordinary port connected to a concatenation
// already splits into one alias segment per operand, so this is the
// same split reached through the explicit spelling.
module explicit_port_concat_child (
    input  logic [7:0] d,
    output logic [7:0] o
);
  assign o = d;
endmodule

module explicit_port_concat (
    input logic clk,
    input logic [7:0] d,
    output .o({x, y})
);

  logic [3:0] x, y;

  explicit_port_concat_child u (.d(d), .o({x, y}));

  assert property (@(posedge clk) x == d[7:4] && y == d[3:0]);

endmodule
