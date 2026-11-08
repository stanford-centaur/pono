// A property whose inductive invariant has to be discovered rather
// than unrolled to.
//
// `cnt` counts 0 to 39 and wraps, so it never reaches 63 and `bad` is
// never set. Nothing about that is visible one step at a time: from
// the unreachable state cnt == 40 there is a chain of twenty-four
// states, every one with `bad` low, ending at cnt == 63 where `bad`
// is set. k-induction has to unroll past that whole chain before its
// simple-path constraint rules it out -- measured, it reports nothing
// at k = 24 and proves at k = 30 -- while IC3 finds `cnt <= 39` and
// is done.
//
// (k-induction with simple paths is complete for a finite system, so
// this is a difference in what it costs rather than in what is
// possible. The chain is what sets the cost, and widening `cnt`
// widens it.)
module unreachable_window (
    input logic clk
);

  logic [5:0] cnt;
  logic       bad;

  initial cnt = 6'd0;
  initial bad = 1'b0;

  always_ff @(posedge clk) begin
    if (cnt == 6'd39) cnt <= 6'd0;
    else cnt <= cnt + 6'd1;
    if (cnt == 6'd63) bad <= 1'b1;
  end

  assert property (@(posedge clk) !bad);

endmodule
