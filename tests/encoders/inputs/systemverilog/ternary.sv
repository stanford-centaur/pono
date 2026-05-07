// Tests the ?: ternary operator.
// next-state of x picks between two values based on sel.
// Property: x != 9.  BMC must falsify within 2 steps with sel=1.
module ternary (input logic clk, input logic sel);
  logic [3:0] x;
  always_ff @(posedge clk) x <= sel ? 9 : 2;
  assert property (@(posedge clk) x != 9);
endmodule
