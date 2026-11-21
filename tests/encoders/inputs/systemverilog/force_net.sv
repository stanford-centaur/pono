// A `force` on a net rather than a variable. The two differ in what
// `release` does -- a net returns to its drivers' value at once,
// while a variable keeps the forced value until its next procedural
// assignment -- so a net is the simpler of the two to model.
//
// Never released here, so `w` follows the force rather than the
// continuous assign.
module force_net (
    input logic clk,
    input logic d
);

  wire w;
  assign w = d;

  initial force w = 1'b1;

  assert property (@(posedge clk) w == 1'b1);

endmodule
