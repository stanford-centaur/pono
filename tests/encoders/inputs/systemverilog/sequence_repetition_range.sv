// A range of counts over a whole sequence. `[*1:2]` admits the
// single match, which subsumes the double one, so the whole thing
// completes exactly when one `a ##1 b` does.
module sequence_repetition_range (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b
);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> (((a ##1 b)[*1:2]) |-> (b && $past(a))));

endmodule
