// Procedural `disable` of the enclosing named block, used the way
// `break` is used in break_in_for.sv but via a completely different
// statement kind.  Same din == 4'b1000 constant and same
// distinguishing logic: correct semantics disable `main_blk` (which
// also aborts the for loop) before bit 3 is visited, so any_set stays
// 0 forever.
module disable_named_block (input logic clk, input logic rst);

  logic [3:0] din;
  logic       any_set, reg_any_set;

  assign din = 4'b1000;

  always_comb begin : main_blk
    any_set = 1'b0;
    for (int i = 0; i < 4; ++i) begin
      if (i == 2) disable main_blk;
      any_set |= din[i];
    end
  end

  always_ff @(posedge clk) begin
    if (rst) reg_any_set <= 1'b0;
    else reg_any_set <= any_set;
  end

  assert property (@(posedge clk) reg_any_set != 1'b1);

endmodule
