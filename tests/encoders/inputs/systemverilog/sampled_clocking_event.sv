// A sampled-value function's optional `clocking_event` argument,
// naming the clock the design already runs on. Every sample here is
// taken on that clock anyway, so saying so is redundant rather than
// wrong: the argument is checked and dropped.
//
// `$past`'s third argument is its gating `enable`, and omitting it
// while supplying the event after it -- the `, ,` below -- has to
// read as "not supplied" rather than as an expression to evaluate.
module sampled_clocking_event (
    input logic clk,
    input logic a
);

  assert property (@(posedge clk)
                   ($rose(a, @(posedge clk)) == $rose(a))
                   and ($fell(a, @(posedge clk)) == $fell(a))
                   and ($stable(a, @(posedge clk)) == $stable(a))
                   and ($changed(a, @(posedge clk)) == $changed(a))
                   and ($past(a, 1, , @(posedge clk)) == $past(a)));

endmodule
