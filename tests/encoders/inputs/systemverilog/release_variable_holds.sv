// What `release` does to a *variable*: it does not restore anything.
// The variable keeps the value the force left it with until its next
// procedural assignment (LRM 10.6.2), and here it never gets one --
// nothing else assigns `x` at all.
//
// Forced at phase 0 and released at phase 1, so by phase 3 the
// release is two cycles past and `x` still carries the forced value.
module release_variable_holds (
    input logic clk
);

  logic [1:0] phase;
  initial phase = 2'd0;
  always_ff @(posedge clk) if (phase != 2'd3) phase <= phase + 2'd1;

  logic [3:0] x;

  always @(posedge clk) begin
    if (phase == 2'd0) force x = 4'd5;
    else if (phase == 2'd1) release x;
  end

  assert property (@(posedge clk) phase == 2'd3 |-> x == 4'd5);

endmodule
