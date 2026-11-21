// Two accumulators with the identical update rule, one declared
// `bit` (2-state) and one `logic` (4-state).  Pono's SMT bitvector
// model has no X/Z states, so the two should be numerically
// indistinguishable at every cycle -- this genuinely holds forever.
// Paired with bit_vs_logic_parity_mismatch.sv, which asserts the
// opposite and is falsified trivially, so the pair together actually
// exercises both directions of the comparison rather than only ever
// proving equality (or only ever proving inequality) in isolation.
module bit_vs_logic_parity (input logic clk, input logic rst);

  bit   [7:0] acc_bit;
  logic [7:0] acc_logic;

  always_ff @(posedge clk) begin
    if (rst) begin
      acc_bit   <= 8'd0;
      acc_logic <= 8'd0;
    end else begin
      acc_bit   <= acc_bit + 8'd7;
      acc_logic <= acc_logic + 8'd7;
    end
  end

  assert property (@(posedge clk) acc_bit == acc_logic);

endmodule
