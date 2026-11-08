// star_trailing.sv's completion condition shifted one cycle late --
// the error an empty match invites, since it ends where it began.
// Must be refuted.
module star_trailing_fails (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b,
    input logic c
);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2)
                                   && !$past(rst, 3))
                   |-> ((b ##1 a[*]) |-> $past(b || ($past(b) && a))));

endmodule
