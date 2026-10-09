// Holds counterpart to eventually_range_fails.sv.  `count` wraps
// 0,1,2,3,0,... so any four consecutive cycles cover every value:
// at every cycle t, count is 3 somewhere in t+1 .. t+4.
module eventually_range (
    input logic clk,
    input logic rst
);

  logic [1:0] count;

  always_ff @(posedge clk) begin
    if (rst) count <= 2'd0;
    else count <= count + 2'd1;
  end

  assert property (@(posedge clk) eventually [1:4] (count == 2'd3));

endmodule
