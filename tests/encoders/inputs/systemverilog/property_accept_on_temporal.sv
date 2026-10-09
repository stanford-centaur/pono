// An abort over a genuinely temporal operand. The true-padded
// truncation forgives the operand's eventualities from the abort
// onward while still holding it to its safety part before -- a split
// of the property this encoder does not have, since the tableau
// gives one term for both. Refused rather than approximated in
// either direction, both of which are wrong answers rather than
// missing ones.
module property_accept_on_temporal (
    input logic clk,
    input logic abrt,
    input logic a
);

  assert property (@(posedge clk) accept_on (abrt) (s_eventually a));

endmodule
