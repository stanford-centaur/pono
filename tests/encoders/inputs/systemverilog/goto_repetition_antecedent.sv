// A goto count as an implication's antecedent. The antecedent has to
// say its match *ends now*, which for occurrences that need not be
// adjacent is a fact about unbounded history rather than the
// eventuality the tableau would read it as -- but unbounded history
// is not unbounded state. Attempts start at every cycle, and as the
// start moves back the number of occurrences in the window grows one
// at a time, so every count from 1 up to the running total is
// achievable: a match of exactly n exists precisely when the total
// has reached n. A saturating counter says that in
// ceil(log2(n + 1)) bits, which keeps the whole implication a plain
// safety property.
//
// `b` here is at cnt 2, 5 and 9, so the 2nd occurrence -- the first
// cycle the antecedent fires -- is at cnt == 5. The consequent is
// false there, so this is refuted, and goto_repetition_antecedent_
// holds.sv pins the cycle it fires on.
module goto_repetition_antecedent (
    input logic clk,
    input logic rst
);

  logic [3:0] cnt;
  logic b;

  always_ff @(posedge clk) cnt <= rst ? 4'd0 : cnt + 4'd1;
  // Gated on reset: cnt is free at cycle 0, so without this the
  // trace could open with an extra occurrence of b.
  assign b = !rst && ((cnt == 4'd2) || (cnt == 4'd5) || (cnt == 4'd9));

  assert property (@(posedge clk) (!rst) |-> (b[->2] |-> 1'b0));

endmodule
