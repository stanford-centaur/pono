// `$past(a)` sits inside a `&&` operand.  The encoder has to convert
// that operand exactly once: because the history chain a $past call
// builds is not memoized, converting the operand twice would silently
// add a second, redundant chain of latches tracking the same value --
// which no property verdict would reveal, hence the latch-count check
// in the accompanying test.  Holds trivially (the consequent is one of
// the antecedent's own conjuncts).
module past_in_logical_and (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic b
);

  assert property (@(posedge clk) (($past(a) == a) && b) |-> b);

endmodule
