// A range of delays after the empty-admitting repetition. Each delay
// in the range spends one of its cycles on the empty match, so the
// shortest one already reduces `a[*] ##1 b` to plain `b` -- which
// subsumes every longer alternative, and the sequence completes
// exactly when `b` holds.
//
// The one-cycle-late reading is refuted by star_trailing_fails.sv
// and star_middle_fails.sv, which share this fixture's shape.
module star_delay_range (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b,
    input logic c
);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2)
                                   && !$past(rst, 3))
                   |-> ((a[*] ##[1:3] b) |-> b));

endmodule
