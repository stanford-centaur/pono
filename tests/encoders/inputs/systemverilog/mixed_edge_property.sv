// Two properties clocked on opposite edges of the *same* signal
// (`@(posedge clk)` and `@(negedge clk)`). Per the documented design
// decision at the top of assertion_walker.cpp, this encoder has no
// clock-domain-crossing model, so "same clock" isn't enough --
// check_clock() also requires every property to use the same edge,
// and rejects this design too.
module mixed_edge_property (
    input logic clk,
    input logic a,
    input logic b
);

  assert property (@(posedge clk) a |-> b);
  assert property (@(negedge clk) b |-> a);

endmodule
