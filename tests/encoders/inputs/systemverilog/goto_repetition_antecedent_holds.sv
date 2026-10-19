// The counterpart to goto_repetition_antecedent.sv: a `[->n]` match
// always lands *on* an occurrence, unlike `[=n]`, which may end at
// any later cycle before the next one. Proving this is what rules
// out a counter that has simply latched high and fires everywhere.
module goto_repetition_antecedent_holds (
    input logic clk,
    input logic rst,
    input logic b
);

  assert property (@(posedge clk) b[->2] |-> b);

endmodule
