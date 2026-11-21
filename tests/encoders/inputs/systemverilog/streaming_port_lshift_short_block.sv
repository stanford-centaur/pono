// The same cross-check where the slice does *not* divide the width
// (12 bits in blocks of 5, leaving 2). LRM 11.4.14.2 cuts blocks
// from the right-hand end, so the short one starts out leftmost and
// the reversal leaves it at the right -- the case an interval
// calculation gets subtly wrong, and the reason the segments are
// built by walking bits and coalescing.
module st_port_src (input logic [11:0] d, output logic [11:0] o);
  assign o = d;
endmodule

module streaming_port_lshift_short_block (
    input logic clk,
    input logic [11:0] d
);

  logic [3:0] a, b, c;
  logic [3:0] a2, b2, c2;

  st_port_src u (.d(d), .o({<<5{a, b, c}}));

  always_comb {<<5{a2, b2, c2}} = d;

  assert property (@(posedge clk) a == a2 && b == b2 && c == c2);

endmodule
