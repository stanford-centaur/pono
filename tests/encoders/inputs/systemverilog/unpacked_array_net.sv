// An unpacked-array net driven by continuous assigns. A net is not a
// bit range, so driving one is a constraint rather than an
// assignment: each driven element equals its driver, and the array
// itself starts as a free variable so anything left undriven reads
// as an unknown -- which is what an undriven net is.
//
// `m` is driven element by element, `d` is the same over a
// descending range (where the index has to be normalized the other
// way), and `w` is driven whole. For a whole-array target the assign
// *is* the wire's definition, since a wire has no term until its
// driver supplies one.
//
// unpacked_array_net_fails.sv pins the undriven half.
module unpacked_array_net (
    input logic clk,
    input logic rst,
    input logic [7:0] v
);

  wire [7:0] m[0:3];
  wire [7:0] d[3:0];
  logic [7:0] src[0:1];
  wire [7:0] w[0:1];

  assign m[0] = v;
  assign m[1] = v + 8'd1;
  assign m[2] = 8'd9;

  assign d[3] = v;
  assign d[0] = v + 8'd1;

  assign src[0] = v;
  assign src[1] = v + 8'd2;
  assign w = src;

  assert property (@(posedge clk)
                   m[0] == v && m[1] == v + 8'd1 && m[2] == 8'd9
                   && d[3] == v && d[0] == v + 8'd1
                   && w[0] == v && w[1] == v + 8'd2);

endmodule
