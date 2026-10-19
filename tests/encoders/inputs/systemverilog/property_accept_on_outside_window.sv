// The negative of property_accept_on_window.sv: the same design with
// the pulse moved clear of the failing attempt, so nothing is waived
// and the failure stands. Without it, an abort that simply waived
// everything would pass the positive.
module property_accept_on_outside_window (input logic clk);

  logic [3:0] cnt;
  initial cnt = 4'd0;
  always @(posedge clk) cnt <= cnt + 4'd1;

  logic abrt, a, b;
  assign abrt = (cnt == 4'd7);
  assign a   = (cnt == 4'd0);
  assign b   = 1'b0;

  assert property (@(posedge clk) accept_on (abrt) (a |=> b));

endmodule
