// Same request/response liveness as req_ack.sv, but here `ack` is tied
// high, so every request is acknowledged in the same cycle and
// G(req -> F ack) genuinely holds.  There is no fair lasso violating
// it, so the LTL tableau must NOT manufacture a counterexample: BMC
// exhausts its bound without finding a witness and returns UNKNOWN.
// This guards against an unsound tableau that reports spurious
// liveness failures.
module req_ack_holds (input logic clk, input logic rst,
                      input logic req);
  logic ack;
  assign ack = 1'b1;
  assert property (@(posedge clk) always (req |-> s_eventually ack));
endmodule
