// A `force` applied once and never released, so the override is in
// effect for the whole run. Nothing here is measured in simulation
// time: the extent is "always".
//
// `force` drives its target continuously and overrides every other
// driver (LRM 10.6.2), so the counter's assignments do not reach `x`.
module force_permanent (
    input logic clk
);

  logic [3:0] x;
  always_ff @(posedge clk) x <= x + 4'd1;

  initial force x = 4'd5;

  assert property (@(posedge clk) x == 4'd5);

endmodule
