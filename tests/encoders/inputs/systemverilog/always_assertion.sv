// `always P` at the top of an SVA assertion is treated as safety
// -- the encoder unwraps it and emits P as a per-cycle invariant,
// exactly as it would for a plain `assert property (@... P)`.
// Counter increments on every cycle; the assertion fails when it
// reaches 4 (one reset cycle + four increments).
module always_assertion (input logic clk, input logic rst);
  logic [4:0] count;
  always_ff @(posedge clk) count <= rst ? '0 : count + 1;
  assert property (@(posedge clk) always (count != 4));
endmodule
