// The bounded twin of unbounded_repeat_plus.sv: plain `a`, refuted
// at the same depth.
module unbounded_repeat_plus_bounded (
    input logic clk,
    input logic rst,
    input logic a
);

  assert property (@(posedge clk) a |-> 1'b0);

endmodule
