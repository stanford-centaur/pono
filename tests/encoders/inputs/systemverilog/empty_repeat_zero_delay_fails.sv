// The other side of empty_repeat_zero_delay.sv: dropping the empty
// branch must not also drop the one-occurrence branch. `a ##0 b` --
// a single `a` overlapping `b` -- is still a match, so requiring a
// second `a` behind it is too strong and must be refuted.
module empty_repeat_zero_delay_fails (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b
);

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> ((a[*0:2] ##0 b) |-> (a && b && $past(a))));

endmodule
