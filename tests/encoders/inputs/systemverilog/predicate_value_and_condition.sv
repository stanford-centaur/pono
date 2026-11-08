// Each comparison appears in both roles the encoder treats
// differently: as a *value* (assigned into a 1-bit net, so the
// predicate has to be materialised as a bit-vector) and as a
// *condition* (an operand of `||`/`==`, where the encoder wants a
// native Bool).  `<=` is exactly `<` or `==`, so the property is a
// tautology that fails only if the two forms of the same comparison
// ever disagree.  Paired with predicate_value_and_condition_fails.sv.
module predicate_value_and_condition (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);

  logic le_val, lt_val, eq_val;

  assign le_val = (a <= b);
  assign lt_val = (a < b);
  assign eq_val = (a == b);

  logic reg_le, reg_split;

  always_ff @(posedge clk) begin
    if (rst) begin
      reg_le    <= 1'b0;
      reg_split <= 1'b0;
    end else begin
      // Left: the comparison's value, straight off the net.
      reg_le    <= le_val;
      // Right: the same comparisons as conditions, recombined.
      reg_split <= (a < b) || (a == b);
    end
  end

  assert property (@(posedge clk)
                    reg_le == reg_split && le_val == (lt_val || eq_val));

endmodule
