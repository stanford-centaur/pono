// The discriminating negative for the three above: `<<` on the port
// against `>>` procedurally. The two orderings differ, so this must
// be refuted -- and it is what catches a port path that quietly
// ignored the slice and split the bits as `>>` does, which every
// positive test here would still pass.
module st_port_src (input logic [11:0] d, output logic [11:0] o);
  assign o = d;
endmodule

module streaming_port_lshift_mismatch (
    input logic clk,
    input logic [11:0] d
);

  logic [3:0] a, b, c;
  logic [3:0] a2, b2, c2;

  st_port_src u (.d(d), .o({<<4{a, b, c}}));

  always_comb {>>4{a2, b2, c2}} = d;

  assert property (@(posedge clk) a == a2 && b == b2 && c == c2);

endmodule
