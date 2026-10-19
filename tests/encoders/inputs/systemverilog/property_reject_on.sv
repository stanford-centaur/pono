// `reject_on` with the abort condition tied high and the operand
// tied true. An accept waives a failure; a reject *causes* one, so
// this must be refuted even though the operand never fails.
module property_reject_on (input logic clk);

  logic abrt, a;
  assign abrt = 1'b1;
  assign a   = 1'b1;

  assert property (@(posedge clk) reject_on (abrt) a);

endmodule
