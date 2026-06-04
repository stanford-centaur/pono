// Element-select LHS.  A 4-bit `flag` register has each bit set
// individually by a procedural for loop -- `flag[i] <= 1'b1;` --
// inside an always_ff block.  Replaces bit `i` of the prior pending
// next-state value with the new slice, so after one reset cycle
// the register holds 4'b1111 and the assertion fails.
module element_select_lhs (input logic clk, input logic rst);
  logic [3:0] flag;

  always_ff @(posedge clk) begin
    if (rst) begin
      flag <= '0;
    end else begin
      for (int i = 0; i < 4; ++i) begin
        flag[i] <= 1'b1;
      end
    end
  end

  // After the reset cycle, flag becomes 4'b1111 in one more cycle.
  // assert flag != 4'b1111 fails at cycle 2.
  assert property (@(posedge clk) flag != 4'b1111);
endmodule
