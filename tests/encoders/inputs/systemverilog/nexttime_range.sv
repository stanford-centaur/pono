// `nexttime [3]` must shift three cycles, not one.  `count` wraps
// 0,1,2,3,0,..., so three cycles after a 0 it is 3.  An encoder that
// drops the [3] and emits a single `X` would check count == 3 one
// cycle after the 0, where it is actually 1, and report a violation.
// Paired with nexttime_range_fails.sv, which asks for the wrong
// distance.
module nexttime_range (
    input logic clk,
    input logic rst
);

  logic [1:0] count;

  always_ff @(posedge clk) begin
    if (rst) count <= 2'd0;
    else count <= count + 2'd1;
  end

  assert property (@(posedge clk)
                    always ((!rst && count == 2'd0)
                            |-> nexttime [3] (count == 2'd3)));

endmodule
