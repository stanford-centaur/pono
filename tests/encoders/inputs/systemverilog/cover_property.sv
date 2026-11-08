// `cover property` alone (no assert/assume in the design).  Unlike
// assume/restrict, dropping cover silently isn't dangerous for a
// safety-focused model checker -- there's no proof obligation to
// corrupt -- but it should still be *cleanly* absent rather than
// miscounted as a safety property or a liveness obligation.
module cover_property (input logic clk, input logic rst, input logic [3:0] data);

  cover property (@(posedge clk) data == 4'd5);

endmodule
