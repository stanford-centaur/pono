// A sampled-value function naming a *different* clock than the one
// the design runs on. The sampled-value path and the property path
// consult the same record of the design's one clock, so both refuse
// it.
module sampled_clocking_event_second_clock (
    input logic clk,
    input logic clk2,
    input logic a
);

  assert property (@(posedge clk) $rose(a, @(posedge clk2)) == $rose(a));

endmodule
