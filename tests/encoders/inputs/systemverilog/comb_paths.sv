// The same XOR function computed three different ways -- a
// continuous `assign`, an `always_comb` block, and a legacy
// `always @*` block -- each latched into its own register so the
// comparison is over consistently-paired same-cycle snapshots.
// Genuinely holds for every a, b, since all three compute the same
// function.  Paired with comb_paths_fails.sv.
module comb_paths (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);

  logic [3:0] via_assign, via_comb, via_legacy;
  logic [3:0] via_assign_reg, via_comb_reg, via_legacy_reg;

  assign via_assign = a ^ b;

  always_comb via_comb = a ^ b;

  always @* via_legacy = a ^ b;

  always_ff @(posedge clk) begin
    if (rst) begin
      via_assign_reg <= 4'd0;
      via_comb_reg   <= 4'd0;
      via_legacy_reg <= 4'd0;
    end else begin
      via_assign_reg <= via_assign;
      via_comb_reg   <= via_comb;
      via_legacy_reg <= via_legacy;
    end
  end

  assert property (@(posedge clk)
                    via_assign_reg == via_comb_reg && via_comb_reg == via_legacy_reg);

endmodule
