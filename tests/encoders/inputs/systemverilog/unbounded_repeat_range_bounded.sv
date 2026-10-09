// The bounded twin of unbounded_repeat_range.sv.
module unbounded_repeat_range_bounded (
    input logic clk,
    input logic rst,
    input logic a
);

  assert property (@(posedge clk) a[*2] |-> 1'b0);

endmodule
