// predicate_value_and_condition.sv with `<` swapped for `>` on the
// condition side, so the two forms of the comparison genuinely differ
// and the property is violated as soon as a > b.  Confirms the holds
// variant is actually checking the value/condition agreement rather
// than passing vacuously.
module predicate_value_and_condition_fails (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);

  logic le_val, gt_val, eq_val;

  assign le_val = (a <= b);
  assign gt_val = (a > b);
  assign eq_val = (a == b);

  assert property (@(posedge clk) le_val == (gt_val || eq_val));

endmodule
