// A procedural *immediate* assertion (`assert (expr);`), distinct
// from the concurrent `assert property (...)` form used everywhere
// else in this suite.  The encoder header's doc-comment claims
// "Basic SVA immediate assertions" are supported -- this is the test
// that actually exercises that claim.  Both `rst` and `a` are plain
// (unregistered) inputs, already free from cycle 1 onward, so BMC can
// pick rst == 0 and a == 7 simultaneously at cycle 1 itself.
module immediate_assert (input logic clk, input logic rst, input logic [3:0] a);

  always_ff @(posedge clk) begin
    if (!rst) begin
      assert (a != 4'd7);
    end
  end

endmodule
