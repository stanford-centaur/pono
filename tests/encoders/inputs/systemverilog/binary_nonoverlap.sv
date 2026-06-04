// SVA non-overlapping implication (`|=>`).  At every cycle, if
// `arm` was high last cycle, the registered count must be
// non-zero this cycle.  Implementation introduces a hidden latch
// state var so "arm was high last cycle" is current-state-only.
//
// arm is a free input.  cnt resets to 0 and otherwise increments
// when arm is high (so it transitions 0->1 in the cycle after arm
// asserts).  With reset preprocessing, the antecedent first fires
// at cycle 1 (arm=1), the consequent fires at cycle 2 -- where
// the implementation's `cnt` is still 0 if arm wasn't asserted
// last cycle.  Falsifiable: arm=0 at cycle 1, arm=1 at cycle 2;
// cycle 3 checks cnt != 0 with cnt=1.  Property holds.  To
// falsify we need arm=1 last cycle AND cnt==0 this cycle, i.e.
// rst=0, arm=1 at cycle 1, then cnt=1 by cycle 2 ... actually
// this property is true by construction.
//
// Use a simpler falsifiable property instead: if arm was high
// last cycle, data this cycle must equal 4'b1010.  data is a
// free input; BMC picks data != 4'b1010 at cycle 2 (after arm=1
// at cycle 1) to falsify.
module binary_nonoverlap (input logic clk, input logic rst,
                          input logic arm, input logic [3:0] data);
  // The consequent fires one cycle after arm.  Reset preprocessing
  // keeps cycle 0 vacuous; arm=1 at cycle 1 sets the latch; at
  // cycle 2 the property checks (data == 4'b1010).  BMC picks
  // data = anything-but-10 at cycle 2 -> falsifies.
  assert property (@(posedge clk) arm |=> data == 4'b1010);
endmodule
