// Integration test #3: a minimal credit-based flow-control link
// (hierarchical sender instance) checked with IC3Bits rather than
// Bmc/KInduction, to exercise a third distinct prover engine against
// the SV frontend's output.  This design exists specifically to show
// the assume_property.sv gap (test_systemverilog_sva.cpp) mattering
// at design scale, not just in a single-property unit test: the
// environment assumption ("never push when out of credits") is
// exactly the kind of backpressure contract real hierarchical designs
// rely on `assume property` to state, and its being a complete no-op
// today means IC3 finds a "violation" that a correctly-behaving
// environment could never actually trigger.
// Deliberately unguarded: without the (ignored) `assume property`
// below, nothing stops `push` while `credits == 0`, so credits
// underflows to 4'hF.  A real sender would rely on the environment
// contract the assume statement is meant to express, rather than
// re-deriving a saturating guard here; this mirrors how the
// assume_property.sv unit test also puts all the "protection" in the
// assumption rather than the design itself.
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
