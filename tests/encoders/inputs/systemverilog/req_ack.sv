// Canonical request/response liveness: every request must eventually
// be acknowledged.  This nests `always` over `s_eventually`
// (G(req -> F ack)), so it is NOT expressible as a pure safety
// property and exercises the general LTL tableau.  `ack` is a free
// input, so BMC can withhold it forever after a request -- a fair
// lasso along which the F-eventuality is never discharged is exactly
// the counterexample.
module req_ack (input logic clk, input logic rst,
                input logic req, input logic ack);
  assert property (@(posedge clk) always (req |-> s_eventually ack));
endmodule
