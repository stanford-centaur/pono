// `s_eventually P`: a liveness property claiming P must hold at
// some future cycle.  The encoder pushes !P onto liveness_propvec_;
// the test wraps the resulting TS with LivenessToSafetyTranslator
// so BMC can hunt for the lasso counterexample.
//
// Here `enable` is a free input.  When `enable` is held low after
// the reset, `count` stays at 0 forever and never reaches 5, so
// the property is violated -- the L2S-translated safety property
// becomes refutable.
module eventually_assertion (input logic clk, input logic rst,
                             input logic enable);
  logic [4:0] count;
  always_ff @(posedge clk) begin
    if (rst)         count <= '0;
    else if (enable) count <= count + 1;
  end
  assert property (@(posedge clk) s_eventually (count == 5));
endmodule
