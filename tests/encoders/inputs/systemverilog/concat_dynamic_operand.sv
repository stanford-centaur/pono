// A concatenation target with one runtime-indexed operand. The
// dynamic operand cannot be routed through the bit-range write path,
// and dropping it would drop `hi`'s write with it, so the whole
// assignment is rejected rather than half-encoded.
module concat_dynamic_operand (
    input logic clk,
    input logic [1:0] i,
    input logic [1:0] src
);

  logic hi;
  logic [3:0] arr;

  always_ff @(posedge clk) {hi, arr[i]} <= src;

  assert property (@(posedge clk) hi == hi);

endmodule
