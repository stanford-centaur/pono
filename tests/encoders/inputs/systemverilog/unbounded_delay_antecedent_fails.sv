// unbounded_delay_antecedent.sv's condition shifted a cycle late.
// A latch that remembered the wrong cycle would let this hold.
module unbounded_delay_antecedent_fails (
    input logic clk,
    input logic rst,
    input logic ai,
    input logic bi
);

  logic a, b, seen_a;

  assign a = ai && !rst;
  assign b = bi && !rst;
  always_ff @(posedge clk) seen_a <= rst ? 1'b0 : (seen_a || a);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> ((a ##[1:$] b) |-> $past(b && seen_a)));

endmodule
