// Tests `initial` blocks: x must start at 7 and never change.
// Property: x != 7.  BMC must falsify at step 0 because the initial
// state has x = 7.
module initial_block (input logic clk);
  logic [3:0] x;
  initial x = 7;
  always_ff @(posedge clk) x <= x;
  assert property (@(posedge clk) x != 7);
endmodule
