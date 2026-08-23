// A range-select lvalue with a non-constant (runtime-variable) base
// (`w[base +: 4]`, the indexed part-select form, which the LRM allows
// with a variable base as long as the width is constant). There is no
// dynamic-range-select write fallback anywhere in this encoder (unlike
// ElementSelect's single-bit dynamic-index fallback), so this must
// throw a clear PonoException rather than silently dropping the write.
module dynamic_range_select_lhs (input logic clk,
                                 input logic [2:0] base,
                                 input logic [3:0] a);

  logic [7:0] w;

  always_ff @(posedge clk) begin
    w[base+:4] <= a;
  end

endmodule
