// The reading unbounded_repeat_star.sv collapses to: the consequent
// on its own. Sharing a refutation depth with it is what shows the
// empty antecedent really is unconditional rather than merely
// encodable.
module unbounded_repeat_star_equiv (
    input logic clk,
    input logic rst
);

  logic [3:0] cnt;
  logic c;

  always_ff @(posedge clk) cnt <= rst ? 4'd0 : cnt + 4'd1;
  assign c = !rst && (cnt != 4'd5);

  assert property (@(posedge clk) (!rst) |-> c);

endmodule
