// empty_repeat_match.sv with its expansion shifted one cycle late.
// The empty match ends where it began, so a reading that puts it a
// cycle later is exactly the off-by-one the LRM's delay-absorption
// rule rules out -- this must be refuted.
module empty_repeat_match_fails (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b
);

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> ((a[*0:2] ##1 b) |-> $past(b)));

endmodule
