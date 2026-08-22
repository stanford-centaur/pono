// Range-select (bit-slice) LHS on a plain continuous assign: two
// separate assigns each pin a different half of `w` via a constant
// range select (`w[3:0]`, `w[7:4]`), rather than a single-bit
// ElementSelect. resolve_lvalue() must recognize RangeSelect the same
// way it already recognizes ElementSelect/MemberAccess so these
// constrain `w` instead of being silently dropped.
module range_select_lhs (input logic clk,
                          input logic [3:0] a,
                          input logic [3:0] b);

  logic [7:0] w;

  assign w[3:0] = a;
  assign w[7:4] = b;

  assert property (@(posedge clk) w == { b, a });

endmodule
