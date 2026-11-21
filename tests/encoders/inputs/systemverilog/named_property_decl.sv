// A named `sequence`/`property` declaration, referenced by name from
// `assert property (p_check)` -- every other assertion in this suite
// uses an inline property expression instead.  s_arm matches exactly
// when data == 5; p_check's (deliberately self-contradictory)
// consequent means the property is falsified the moment data == 5,
// which `data` (free) can pick at cycle 1.
module named_property_decl (input logic clk, input logic rst, input logic [3:0] data);

  sequence s_arm;
    @(posedge clk) data == 4'd5;
  endsequence

  property p_check;
    @(posedge clk) s_arm |-> data != 4'd5;
  endproperty

  assert property (p_check);

endmodule
