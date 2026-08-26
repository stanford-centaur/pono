// Integration test #1 (liveness half): identical hierarchy/arbitration
// logic to arbiter_bus.sv, but with a liveness property in place of the
// safety one, checked via the LivenessToSafetyTranslator path
// (check_liveness_bmc) instead of plain Bmc.
//
// "Whenever requester 0 asks, it must eventually be granted."  `req`
// is free, so BMC/L2S can pick req[1] held continuously while req[0]
// pulses once, starving requester 0 if `turn` never favors it.
typedef struct packed {
  logic valid;
  logic id;
} grant_t;

module requester (input logic clk, input logic rst, input logic req_in, input logic grant_in);
  logic pending;
  always_ff @(posedge clk) begin
    if (rst) pending <= 1'b0;
    else if (req_in) pending <= 1'b1;
    else if (grant_in) pending <= 1'b0;
  end
endmodule

module arbiter_bus_live (input logic clk, input logic rst, input logic [1:0] req);

  logic   turn;
  grant_t grant;

  always_comb begin
    if (req[0] && (!req[1] || !turn)) grant = '{valid: 1'b1, id: 1'b0};
    else if (req[1]) grant = '{valid: 1'b1, id: 1'b1};
    else grant = '{valid: 1'b0, id: 1'b0};
  end

  requester r0 (
      .clk     (clk),
      .rst     (rst),
      .req_in  (req[0]),
      .grant_in(grant.valid && !grant.id)
  );

  requester r1 (
      .clk     (clk),
      .rst     (rst),
      .req_in  (req[1]),
      .grant_in(grant.valid && grant.id)
  );

  always_ff @(posedge clk) begin
    if (rst) turn <= 1'b0;
    else if (grant.valid) turn <= ~grant.id;
  end

  assert property (@(posedge clk)
                    always (req[0] |-> s_eventually (grant.valid && !grant.id)));

endmodule
