// Companion to comb_paths.sv: the legacy `always @*` path computes a
// *different* function (AND instead of XOR) from the other two, so
// the three-way equivalence breaks as soon as free a, b make XOR and
// AND disagree -- the earliest opportunity is cycle 2 (cycle 1 is
// still the all-zero reset value, where XOR and AND happen to agree).
module comb_paths_fails (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);

  logic [3:0] via_comb, via_legacy;
  logic [3:0] via_comb_reg, via_legacy_reg;

  always_comb via_comb = a ^ b;

  always @* via_legacy = a & b;

  always_ff @(posedge clk) begin
    if (rst) begin
      via_comb_reg   <= 4'd0;
      via_legacy_reg <= 4'd0;
    end else begin
      via_comb_reg   <= via_comb;
      via_legacy_reg <= via_legacy;
    end
  end

  assert property (@(posedge clk) via_comb_reg == via_legacy_reg);

endmodule
