// Operands of unequal width against a slice that divides neither of
// them, so operand boundaries and block boundaries fall in different
// places and a single operand is served by several segments.
module st_port_src (input logic [11:0] d, output logic [11:0] o);
  assign o = d;
endmodule

module streaming_port_lshift_uneven (
    input logic clk,
    input logic [11:0] d
);

  logic [6:0] a;
  logic [4:0] b;
  logic [6:0] a2;
  logic [4:0] b2;

  st_port_src u (.d(d), .o({<<3{a, b}}));

  always_comb {<<3{a2, b2}} = d;

  assert property (@(posedge clk) a == a2 && b == b2);

endmodule
