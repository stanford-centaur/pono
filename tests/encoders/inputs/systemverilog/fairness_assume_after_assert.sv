// The assumption written after the assertion it constrains. A
// property is sealed when it is walked, so the fairness conditions
// are distributed in a pass over the whole design rather than as
// each assertion is met; this is the design that tells the two
// apart.
module fairness_assume_after_assert (input logic clk, input logic a);

  assert property (@(posedge clk) s_eventually a);
  assume property (@(posedge clk) s_eventually a);

endmodule
