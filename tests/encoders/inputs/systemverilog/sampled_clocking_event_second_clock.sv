// A sampled-value function naming a *different* clock than the one
// the design runs on. This used to be accepted and ignored, so
// `$rose(a, @(posedge clk2)) == $rose(a)` proved -- the second clock
// read as the first -- while a *property* clocked on that same clk2
// was rejected outright. The two paths now consult the same record of
// the design's one clock, so both refuse it.
module sampled_clocking_event_second_clock (
    input logic clk,
    input logic clk2,
    input logic a
);

  assert property (@(posedge clk) $rose(a, @(posedge clk2)) == $rose(a));

endmodule
