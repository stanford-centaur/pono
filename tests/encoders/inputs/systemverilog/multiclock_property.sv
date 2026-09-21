// A multiclock property (`@(posedge clk1) a ##1 @(posedge clk2) b`).
// Per the documented design decision at the top of
// assertion_walker.cpp, this encoder has no clock-domain-crossing
// model, so a property clocked on more than one signal (or more than
// one edge of the same signal) is rejected outright.
module multiclock_property (
    input logic clk1,
    input logic clk2,
    input logic rst,
    input logic a,
    input logic b
);

  assert property (@(posedge clk1) a ##1 @(posedge clk2) b |-> 1'b0);

endmodule
