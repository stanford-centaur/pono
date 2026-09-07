// `strong(seq)`: turns a sequence into a liveness obligation that it
// must eventually match.  `a` is a free input that BMC can hold low
// forever, so the sequence (which needs `a` true to even start)
// never matches -- a fair lasso along which the obligation is never
// discharged.  Exact bound verified against the actual run.
module strong_seq (input logic clk, input logic rst, input logic a);

  assert property (@(posedge clk) strong(a ##1 1'b1));

endmodule
