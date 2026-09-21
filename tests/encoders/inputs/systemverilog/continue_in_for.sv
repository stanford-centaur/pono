// `continue` inside the already-supported `for` loop.  `din` is a
// fixed constant (4'b0010, only bit 1 set) so the outcome
// unambiguously distinguishes "continue honored" (bit 1's iteration
// is skipped, any_set stays 0, property holds forever) from
// "continue ignored" (any_set goes high, falsifying at cycle 2).
module continue_in_for (input logic clk, input logic rst);

  logic [3:0] din;
  logic       any_set, reg_any_set;

  assign din = 4'b0010;

  always_comb begin
    any_set = 1'b0;
    for (int i = 0; i < 4; ++i) begin
      if (i == 1) continue;
      any_set |= din[i];
    end
  end

  always_ff @(posedge clk) begin
    if (rst) reg_any_set <= 1'b0;
    else reg_any_set <= any_set;
  end

  assert property (@(posedge clk) reg_any_set != 1'b1);

endmodule
