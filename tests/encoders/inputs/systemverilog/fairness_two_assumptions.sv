// Two fairness assumptions, and an assertion that neither grants on
// its own -- both sets of conditions have to reach the property.
module fairness_two_assumptions (
    input logic clk,
    input logic a,
    input logic b
);

  assume property (@(posedge clk) s_eventually a);
  assume property (@(posedge clk) s_eventually b);
  assert property (@(posedge clk) s_eventually (a || b));

endmodule
