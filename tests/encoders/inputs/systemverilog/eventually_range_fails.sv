// `count` wraps 0,1,2,3,0,... so it reaches 3 every fourth cycle.
// `eventually [2:3] (count == 3)` asks, at EVERY cycle t, for count to
// be 3 at t+2 or t+3 -- false at t == 2 (count is 0 then 1 there), so
// the property genuinely fails.  Dropping the [2:3] window turns this
// into plain `s_eventually (count == 3)`, which is TRUE, so an encoder
// that ignores the range reports no violation.  This is the regression
// test for that unsoundness.
module eventually_range_fails (
    input logic clk,
    input logic rst
);

  logic [1:0] count;

  always_ff @(posedge clk) begin
    if (rst) count <= 2'd0;
    else count <= count + 2'd1;
  end

  // The explicit `always` is required: a temporal property reaching
  // the LTL tableau is only obliged at cycle 0, and "at every cycle"
  // is what makes the window bite.
  assert property (@(posedge clk)
                    always (eventually [2:3] (count == 2'd3)));

endmodule
