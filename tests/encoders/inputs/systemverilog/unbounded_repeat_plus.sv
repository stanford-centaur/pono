// `a[+]` is `a[*1:$]`. Ending now, a run of at least one ends with
// a run of exactly one, so as an antecedent it says what plain `a`
// says -- and is refuted at the same depth, which is what
// SequenceRepetitionPlus and its bounded twin below check.
module unbounded_repeat_plus (
    input logic clk,
    input logic rst,
    input logic a
);

  assert property (@(posedge clk) a[+] |-> 1'b0);

endmodule
