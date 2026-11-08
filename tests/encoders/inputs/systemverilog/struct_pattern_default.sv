// The `'{default: v}` form of a structure assignment pattern (LRM
// 10.9.2), which fills every member without naming them. The by-name
// spelling of the same value holds, and so must this.
typedef struct packed {
  logic [3:0] hi;
  logic [3:0] lo;
} sp_default_t;

module struct_pattern_default (input logic clk);

  sp_default_t s;
  always_comb s = '{default: 4'h4};

  assert property (@(posedge clk) s.hi == 4'h4 && s.lo == 4'h4);

endmodule
