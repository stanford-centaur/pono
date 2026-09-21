// An empty match meeting a zero delay. The LRM (16.9.2.1, and
// F.3.4.2.2 for the formal version) gives `empty ##0 seq` and
// `seq ##0 empty` no match at all -- unlike `##1`, where the empty
// match absorbs the delay and the concatenation reduces to the
// other operand. So the empty alternative of `a[*0:2]` drops out
// entirely across a `##0`:
//
//   a[*0:2] ##0 b  ==  (a ##0 b) or (a[*2] ##0 b)  ==  a && b
//   b ##0 a[*0:2]  ==  (b ##0 a) or (b ##0 a[*2])
//
// Reading the empty branch as an identity instead would let plain
// `b` match, which is what empty_repeat_zero_delay_fails.sv rules
// out. Contrast empty_repeat_match.sv, which pins the `##1` rule.
module empty_repeat_zero_delay (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b
);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> (((a[*0:2] ##0 b) |-> (a && b))
                        and ((b ##0 a[*0:2])
                             |-> ((a && b) || (a && $past(a && b))))));

endmodule
