// SVA non-overlapping implication (`|=>`): if `arm` held last cycle,
// `data` this cycle must equal 4'b1010.  Implemented via a 1-cycle
// latch chain delaying the antecedent so the resulting implication is
// current-cycle (make_history_chain() in sva.cpp).
//
// `arm`/`data` are free inputs; BMC sets arm=1 at cycle 1 and
// data != 4'b1010 at cycle 2 to falsify.
module binary_nonoverlap (input logic clk, input logic rst,
                          input logic arm, input logic [3:0] data);
  assert property (@(posedge clk) arm |=> data == 4'b1010);
endmodule
