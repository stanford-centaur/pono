// The by-position form of a structure assignment pattern (LRM
// 10.9.2's first example, `'{1, 2+k}`), which is a different
// expression class from the by-name form. The by-name spelling of
// the same value holds, and so must this.
typedef struct packed {
  logic [3:0] hi;
  logic [3:0] lo;
} sp_pos_t;

module struct_pattern_positional (input logic clk);

  sp_pos_t s;
  always_comb s = '{4'h4, 4'h5};

  assert property (@(posedge clk) s.hi == 4'h4 && s.lo == 4'h5);

endmodule
