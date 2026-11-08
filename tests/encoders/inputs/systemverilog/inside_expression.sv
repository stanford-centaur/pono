// `x inside {...}` -- set membership (LRM 11.4.13). The result is
// the or-reduction of the members' comparisons, where a plain
// member is compared with *wildcard* equality (an x or z bit in the
// set value is a don't-care; one on the left is not) and a
// `[lo:hi]` member is an inclusive range instead.
//
// Each conjunct states the equivalence with the hand-written form,
// so all of them must hold. `4'b1?01` pins bits 3, 1 and 0 and
// leaves bit 2 free, which is why it is `(s & 4'b1011) == 4'b1001`
// rather than an equality. The signed range is there because the
// bounds decide whether the comparison is signed.
module inside_expression (
    input logic clk,
    input logic [3:0] s,
    input logic signed [3:0] t
);

  assert property (@(posedge clk)
                   ((s inside {[4'd1:4'd3]}) == (s >= 4'd1 && s <= 4'd3))
                   && ((s inside {4'd1, 4'd2}) == (s == 4'd1 || s == 4'd2))
                   && ((s inside {4'd0, [4'd2:4'd3]})
                       == (s == 4'd0 || (s >= 4'd2 && s <= 4'd3)))
                   && ((s inside {4'b1?01}) == ((s & 4'b1011) == 4'b1001))
                   && ((s inside {[4'd5:4'd5]}) == (s == 4'd5))
                   && ((t inside {[-4'sd2:4'sd1]})
                       == (t >= -4'sd2 && t <= 4'sd1)));

endmodule
