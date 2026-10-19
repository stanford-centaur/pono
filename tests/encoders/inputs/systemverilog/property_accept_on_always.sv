// The same property with the abort condition tied high: every
// attempt is waived, so it holds however the operand behaves.
// Paired with property_reject_on.sv, which is the same design under
// the opposite action and must be refuted -- the pair is what
// separates the two actions rather than merely exercising both.
module property_accept_on_always (input logic clk, input logic a);

  logic abrt;
  assign abrt = 1'b1;

  assert property (@(posedge clk) accept_on (abrt) a);

endmodule
