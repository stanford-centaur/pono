// A bounded `s_always` window nested under an implication -- the
// position where the window is relative to the match point rather than
// to the start of the trace, so it exercises the forward-shift
// encoding rather than any trace-origin reindexing.
//
// `count` wraps 0,1,2,3,0,..., so the three cycles after a 0 are
// 1,2,3: all nonzero.  Paired with s_always_range_fails.sv, which only
// widens the window by one and so reaches the next 0.
//
// The `!rst` in the antecedent is required: the initial value of
// `count` is free and reset drives it to 0 at the next cycle too, so
// without it the antecedent can match during reset and then see a
// second 0 immediately after -- a real violation, but of the reset
// sequence rather than of what this fixture means to check.
module s_always_range (
    input logic clk,
    input logic rst
);

  logic [1:0] count;

  always_ff @(posedge clk) begin
    if (rst) count <= 2'd0;
    else count <= count + 2'd1;
  end

  // The explicit `always` is required: a temporal property reaching
  // the LTL tableau is only obliged at cycle 0, so without it the
  // reset cycle makes the implication vacuous and nothing is checked.
  assert property (@(posedge clk)
                    always ((!rst && count == 2'd0)
                            |-> s_always [1:3] (count != 2'd0)));

endmodule
