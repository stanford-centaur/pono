// `(seq)[*n]` -- consecutive repetition of a whole sequence, not just
// a Boolean. The LRM makes it `seq` concatenated with itself at
// `##1`, so its offsets are `seq`'s own convolved with themselves n
// times: `(a ##1 b)[*2]` completes exactly when b, a, b, a hold on
// the last four cycles, most recent first.
//
// Stated as a subset claim -- every match the encoder reports is one
// the hand expansion reports -- and paired with
// sequence_repetition_fails.sv, which shifts that expansion a cycle
// so this cannot hold vacuously.
module sequence_repetition (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b
);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2)
                                   && !$past(rst, 3) && !$past(rst, 4))
                   |-> (((a ##1 b)[*2])
                        |-> (b && $past(a) && $past(b, 2) && $past(a, 3))));

endmodule
