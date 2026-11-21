// `a[*]` as the *last* element of a concatenation. `(b ##1 empty)`
// is `b ##0 1`, i.e. `b` on its own, so the empty match pulls the
// end of the sequence back to where `b` matched rather than a cycle
// later. The only other offset is a run of `a` ending now: a longer
// run ending now also ends with exactly one, so `[*]` contributes no
// further alternatives.
//
// Paired with star_trailing_fails.sv.
module star_trailing (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b,
    input logic c
);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2)
                                   && !$past(rst, 3))
                   |-> ((b ##1 a[*]) |-> (b || ($past(b) && a))));

endmodule
