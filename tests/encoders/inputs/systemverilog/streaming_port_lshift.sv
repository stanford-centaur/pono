// A `<<` streaming concatenation as an output-port connection. The
// re-ordering moves bits across the boundaries between the stream's
// expressions, so an operand is no longer one contiguous range of
// the port -- but the blocks it is cut into each are, so the answer
// is more alias segments rather than a different representation.
//
// The check is against the procedural `<<`, which unpacks the same
// value by a separate route (stream_unreorder() on the right-hand
// side, rather than a permutation of alias segments). The two have
// to agree bit for bit; here the slice divides the width, so every
// block is whole.
module st_port_src (input logic [11:0] d, output logic [11:0] o);
  assign o = d;
endmodule

module streaming_port_lshift (input logic clk, input logic [11:0] d);

  logic [3:0] a, b, c;
  logic [3:0] a2, b2, c2;

  st_port_src u (.d(d), .o({<<4{a, b, c}}));

  always_comb {<<4{a2, b2, c2}} = d;

  assert property (@(posedge clk) a == a2 && b == b2 && c == c2);

endmodule
