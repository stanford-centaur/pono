// A register pair written ONLY through a concatenation-target
// nonblocking assignment, with no other write path (no reset branch,
// no plain non-concat assignment anywhere else) -- the pre-scan pass
// that classifies state vars has no Concatenation case, so `carry`/
// `sum` would otherwise never be marked as state vars even though the
// write-processing code fully supports concat-target writes. The
// `initial` block pins a deterministic starting state (mirroring
// initial_block.sv) without going through the buggy classification
// path itself, since plain (non-forever) `initial` blocks are excluded
// from pre-scan entirely.
module concat_lhs_only_write (input logic clk,
                              input logic [3:0] a,
                              input logic [3:0] b);

  logic carry;
  logic [3:0] sum;

  initial begin
    carry = 1'b0;
    sum = 4'd0;
  end

  always_ff @(posedge clk) begin
    { carry, sum } <= { 1'b0, a } + { 1'b0, b };
  end

  // {carry, sum} = {1'b1, 4'd10} = 5'd26 (e.g. a = b = 13) is reachable
  // one cycle after the pinned initial state, once {carry, sum}
  // tracks free inputs a/b via the concat-target write.
  assert property (@(posedge clk) { carry, sum } != 5'd26);

endmodule
