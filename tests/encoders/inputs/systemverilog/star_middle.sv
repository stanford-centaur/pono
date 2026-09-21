// `a[*]` between two other elements, so its empty match has to
// compose on both sides at once: dropping it leaves `b ##1 a ##1 c`,
// and taking it collapses the pair of delays to a single one,
// leaving `b ##1 c`.
//
// Paired with star_middle_fails.sv.
module star_middle (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b,
    input logic c
);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2)
                                   && !$past(rst, 3))
                   |-> ((b ##1 a[*] ##1 c) |-> (c && ($past(b) || ($past(b, 2) && $past(a))))));

endmodule
