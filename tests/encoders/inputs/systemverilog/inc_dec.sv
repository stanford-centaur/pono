// The postfix increment statement (`cnt++;`) as the sole next-state
// update for a register -- distinct from the `cnt <= cnt + 1`
// spelling used elsewhere in this suite, and from the same operator
// used as a `for`-loop step (handled by an entirely separate code
// path). cnt(k) = k-1 for k >= 1, so cnt first equals 3 at cycle 4.
module inc_dec (input logic clk, input logic rst);

  logic [3:0] cnt;

  always_ff @(posedge clk) begin
    if (rst) cnt <= 4'd0;
    else cnt++;
  end

  assert property (@(posedge clk) cnt != 4'd3);

endmodule
