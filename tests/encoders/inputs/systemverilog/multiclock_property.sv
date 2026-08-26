// A multiclock property (`@(posedge clk1) a ##1 @(posedge clk2) b`).
// Per the documented design decision at the top of sva.cpp, this
// encoder has no clock-domain-crossing model: every named clock is
// treated as the same global pono-cycle, so the nested `@(posedge
// clk2)` wrapper is simply unwrapped like any other Clocking node.
module multiclock_property (
    input logic clk1,
    input logic clk2,
    input logic rst,
    input logic a,
    input logic b
);

  assert property (@(posedge clk1) a ##1 @(posedge clk2) b |-> 1'b0);

endmodule
