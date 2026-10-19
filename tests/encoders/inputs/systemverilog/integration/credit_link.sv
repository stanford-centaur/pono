// Integration test #3: a minimal credit-based flow-control link
// (hierarchical sender instance) checked with IC3Bits rather than
// Bmc/KInduction, to exercise a third distinct prover engine against
// the SV frontend's output.
//
// The sender itself is deliberately unguarded: without the
// `assume property` below, nothing stops `push` while `credits == 0`,
// underflowing credits to 4'hF.  The environment assumption ("never
// push when out of credits") is the backpressure contract a real
// sender would rely on instead of re-deriving a saturating guard.
module sender (
    input  logic       clk,
    input  logic       rst,
    input  logic       push,
    output logic [3:0] credits
);
  always_ff @(posedge clk) begin
    if (rst) credits <= 4'd4;
    else if (push) credits <= credits - 4'd1;
  end
endmodule

module credit_link (input logic clk, input logic rst, input logic push);

  logic [3:0] credits;

  sender s (
      .clk    (clk),
      .rst    (rst),
      .push   (push),
      .credits(credits)
  );

  // Environment assumption: never push with zero credits available.
  assume property (@(posedge clk) !(push && credits == 4'd0));

  assert property (@(posedge clk) credits <= 4'd4);

endmodule
