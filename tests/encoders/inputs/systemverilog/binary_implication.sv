// SVA overlapping implication (`|->`).  The property is
// `!rst |-> data != 4'b1111` -- if reset isn't asserted, the
// (free input) data must not be all-ones.  After the one-cycle
// reset, BMC sets data = 4'b1111 to falsify at cycle 1.
module binary_implication (input logic clk, input logic rst,
                           input logic [3:0] data);
  assert property (@(posedge clk) !rst |-> data != 4'b1111);
endmodule
