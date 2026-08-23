// A constant element-select lvalue whose index is out of range for
// its base (`flag[10]` into a 4-bit `flag`). The LRM permits an
// out-of-range constant index (writes are a no-op, reads return 'x),
// but this encoder has no such semantics -- resolve_lvalue() must
// throw a clear PonoException rather than silently computing a
// bit-range that doesn't fit inside `flag`.
module element_select_out_of_bounds_lhs (input logic clk,
                                         input logic a);

  logic [3:0] flag;

  always_ff @(posedge clk) begin
    flag[10] <= a;
  end

endmodule
