// `sync_accept_on`. The synchronous forms differ from the
// asynchronous ones only in being evaluated when the clocking event
// happens rather than at any simulation time step (LRM 24857-24865),
// and this encoder has no time between ticks -- the same argument
// that merges Eventually with SEventually. So this is the same
// property as property_accept_on_window.sv and must give the same
// answer.
module property_sync_accept_on (input logic clk);

  logic [3:0] cnt;
  initial cnt = 4'd0;
  always @(posedge clk) cnt <= cnt + 4'd1;

  logic abrt, a, b;
  assign abrt = (cnt == 4'd0);
  assign a   = (cnt == 4'd0);
  assign b   = 1'b0;

  assert property (@(posedge clk) sync_accept_on (abrt) (a |=> b));

endmodule
