// An antecedent has to say its match ends *now*. For a count of
// occurrences that need not be adjacent that is a fact about
// unbounded history, not an eventuality, so the reading the tableau
// would give it is the wrong one and it is refused instead.
module goto_repetition_antecedent (
    input logic clk,
    input logic a,
    input logic c
);

  assert property (@(posedge clk) a[->2] |-> c);

endmodule
