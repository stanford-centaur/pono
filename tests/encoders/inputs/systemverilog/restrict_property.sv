// `restrict property` is, like `assume property`, meant to constrain
// the state space explored by a formal tool (as opposed to `assume`'s
// broader simulation+formal semantics) -- same test shape as
// assume_property.sv, checking the other constraining-property
// keyword.
module restrict_property (input logic clk, input logic rst, input logic [3:0] data);

  restrict property (@(posedge clk) data == 4'd0);
  assert property (@(posedge clk) data != 4'd7);

endmodule
