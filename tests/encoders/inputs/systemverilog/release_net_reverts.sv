// What `release` does to a *net*, which is the opposite of the
// variable rule: the net is assigned its drivers' value at once (LRM
// 10.6.2), so the forced value does not outlive the release by even
// a cycle.
//
// The contrast with release_variable_holds.sv is the point -- same
// timing, same phases, opposite outcome.
module release_net_reverts (
    input logic clk,
    input logic d
);

  logic [1:0] phase;
  initial phase = 2'd0;
  always_ff @(posedge clk) if (phase != 2'd3) phase <= phase + 2'd1;

  wire w;
  assign w = d;

  always @(posedge clk) begin
    if (phase == 2'd0) force w = 1'b1;
    else if (phase == 2'd1) release w;
  end

  assert property (@(posedge clk) phase == 2'd3 |-> w == d);

endmodule
