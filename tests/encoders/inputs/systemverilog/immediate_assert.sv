// A procedural *immediate* assertion (`assert (expr);`), distinct
// from the concurrent `assert property (...)` form used everywhere
// else in this suite.  The encoder header's doc-comment claims
// "Basic SVA immediate assertions" are supported -- this is the test
// that actually exercises that claim.  `a` is free, so BMC can pick
// a == 7 at the first post-reset transition, falsifying at cycle 2.
module immediate_assert (input logic clk, input logic rst, input logic [3:0] a);

  always_ff @(posedge clk) begin
    if (!rst) begin
      assert (a != 4'd7);
    end
  end

endmodule
