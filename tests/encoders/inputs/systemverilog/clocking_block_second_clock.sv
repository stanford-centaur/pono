// A named clocking block declared on a different clock than the
// design's. Resolving `@(cb)` recurses into the block's own event,
// so the multiclock check sees `clk2` and rejects it exactly as a
// written-out `@(posedge clk2)` would be -- the indirection must
// not become a way around it.
module clocking_block_second_clock (input logic clk, input logic clk2,
                             input logic d);
  logic q;
  clocking cb @(posedge clk2);
  endclocking
  always_ff @(posedge clk) q <= d;
  assert property (@(posedge clk) q == q);
  assert property (@(cb) q == q);
endmodule
