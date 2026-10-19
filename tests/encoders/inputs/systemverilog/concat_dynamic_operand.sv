// A concatenation target one of whose operands is runtime-indexed:
// `{hi, arr[i]} <= src`. The same write on its own, `arr[i] <= x`,
// splices the selected bit onto the variable's own term, so the
// shape has a representation; what it needs is for a concatenation
// operand to reach it.
//
// `src`'s high bit goes to `hi` and its low bit to whichever bit of
// `arr` the index named.
module concat_dynamic_operand (
    input logic clk,
    input logic [1:0] i,
    input logic [1:0] src
);

  logic hi;
  logic [3:0] arr;
  logic started;

  initial started = 1'b0;

  always_ff @(posedge clk) begin
    started <= 1'b1;
    {hi, arr[i]} <= src;
  end

  assert property (@(posedge clk) started
                   |-> (hi == $past(src[1])
                        && arr[$past(i)] == $past(src[0])));

endmodule
