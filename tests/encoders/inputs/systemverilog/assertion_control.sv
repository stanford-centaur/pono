// The `$assert*` family switches assertion checking on and off at
// run time (LRM 20.12). It writes no design state, and still decides
// something this encoder reports on: an assertion the design turned
// off is one whose counterexamples nobody asked for.
module assertion_control (
    input logic clk,
    input logic rst
);

  logic [3:0] cnt;

  initial cnt = 4'd0;

  always_ff @(posedge clk) begin
    if (rst) cnt <= 4'd0;
    else begin
      $assertoff(0);
      cnt <= cnt + 4'd1;
    end
  end

  assert property (@(posedge clk) cnt == 4'd0);

endmodule
