// An empty-admitting repetition inside an `intersect` operand.
// `intersect` matches only where both operands span the same number
// of cycles, and `c ##1 b` spans exactly two, so the empty match --
// which would make the left operand span one -- cannot be used here.
// That the left operand still reports it is what keeps this from
// being over-permissive.
module star_intersect (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b,
    input logic c
);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2)
                                   && !$past(rst, 3))
                   |-> (((a[*] ##1 b) intersect (c ##1 b)) |-> ($past(a) && $past(c) && b)));

endmodule
