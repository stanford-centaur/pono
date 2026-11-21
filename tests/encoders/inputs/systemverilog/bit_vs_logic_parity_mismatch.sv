// Companion to bit_vs_logic_parity.sv: same two identically-updated
// `bit`/`logic` accumulators, but asserting they *differ*.  Both
// start at 0 after reset, so this is falsified at the very first
// post-reset cycle -- the "fails" half of the pair.
module bit_vs_logic_parity_mismatch (input logic clk, input logic rst);

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

  assert property (@(posedge clk) acc_bit != acc_logic);

endmodule
