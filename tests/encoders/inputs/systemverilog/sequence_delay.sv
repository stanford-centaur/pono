// `|-> ##N Q`: at every cycle t, if antecedent holds at t, the
// consequent must hold at t+N.  The encoder builds a chain of N
// 1-bit latches to delay the antecedent so the resulting
// implication is current-cycle.
//
// Property: `arm |-> ##2 data == 4'b1010`.  `arm` is a free
// input; once arm fires at cycle k, two cycles later `data` is
// required to equal 10 -- but `data` is also free, so BMC just
// picks data != 10 to falsify.
//
// Timing:
//   cycle 0: rst=1 (reset), arm free.  Latches 0 by init.
//   cycle 1: rst=0, BMC picks arm=1 -> latch_0 := 1.
//   cycle 2: latch_1 := latch_0 = 1.  data is free.
//   cycle 3: latch_1 == 1; BMC picks data != 10 -> falsify.
//
// So the property fails at cycle 3 (counting reset cycle 0).
module sequence_delay (input logic clk, input logic rst,
                       input logic arm, input logic [3:0] data);
  assert property (@(posedge clk) arm |-> ##2 data == 4'b1010);
endmodule
