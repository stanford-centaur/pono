// Fields must not overlap: writing `a` must leave `b` alone. An
// offset scheme that ran them together would make this hold.
typedef struct {
  logic [3:0] a;
  logic [3:0] b;
} su2_t;

module unpacked_struct_fails (
    input logic clk,
    input logic rst,
    input logic [3:0] d
);

  su2_t s;

  always_ff @(posedge clk) s.a <= d;

  assert property (@(posedge clk) (!rst && !$past(rst)) |-> s.b == $past(d));

endmodule
