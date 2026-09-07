// Strong until as a liveness property.  `busy` is tied high, so
// `busy s_until done` reduces to "done must eventually hold"
// (strong until requires its right operand to occur).  With `done` a
// free input held low forever, the obligation is violated by a lasso
// in which `done` never arrives -- the negated property is the
// release G(!done), enforced from cycle 0.
module until_live (input logic clk, input logic rst, input logic done);
  logic busy;
  assign busy = 1'b1;
  assert property (@(posedge clk) busy s_until done);
endmodule
