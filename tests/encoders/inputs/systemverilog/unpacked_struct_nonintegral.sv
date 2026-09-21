// A struct whose field is itself an unpacked array has no flat
// layout to give that field a bit range in, so it stays refused --
// with a message naming the field rather than the type kind.
typedef struct {
  logic [3:0] a;
  logic [3:0] v[0:1];
} sn_t;

module unpacked_struct_nonintegral (
    input logic clk,
    input logic rst,
    input logic [3:0] d
);

  sn_t s;

  always_ff @(posedge clk) s.a <= d;

  assert property (@(posedge clk) s.a == s.a);

endmodule
