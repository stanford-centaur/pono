// nexttime_range.sv with the wrong shift distance: two cycles after a
// 0 the wrapping counter is 2, not 3, so the property fails.
module nexttime_range_fails (
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
                            |-> nexttime [2] (count == 2'd3)));

endmodule
