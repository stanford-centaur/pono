// `first_match` wrapped directly around the empty-admitting
// repetition, rather than around a concatenation containing one --
// the shape that makes the empty match propagate out through
// `first_match` to the enclosing concatenation. `first_match` only
// restricts which match is reported, never whether one exists, so
// this completes exactly when `b` holds.
module star_first_match (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b,
    input logic c
);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2)
                                   && !$past(rst, 3))
                   |-> ((first_match(a[*]) ##1 b) |-> b));

endmodule
