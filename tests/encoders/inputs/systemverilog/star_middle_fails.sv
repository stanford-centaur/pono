// star_middle.sv's completion condition shifted one cycle late.
// Must be refuted.
module star_middle_fails (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b,
    input logic c
);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2)
                                   && !$past(rst, 3))
                   |-> ((b ##1 a[*] ##1 c) |-> $past(c && ($past(b) || ($past(b, 2) && $past(a))))));

endmodule
