// Concatenation is the inverse of range-select ({hi,lo}[7:4]==hi,
// [3:0]==lo), and a ternary mux always returns whichever operand its
// condition names, all registered from the same cycle's free
// hi/lo/cond so the comparisons are always fairly paired.  Genuinely
// holds for every hi, lo, cond.  Paired with
// ternary_select_concat_fails.sv.
module ternary_select_concat (
    input logic clk,
    input logic rst,
    input logic [3:0] hi,
    input logic [3:0] lo,
    input logic cond
);

  logic [7:0] wide;
  logic [3:0] sel_out, hi_reg, lo_reg;
  logic       cond_reg;

  always_ff @(posedge clk) begin
    if (rst) begin
      wide     <= 8'd0;
      sel_out  <= 4'd0;
      hi_reg   <= 4'd0;
      lo_reg   <= 4'd0;
      cond_reg <= 1'b0;
    end else begin
      wide     <= {hi, lo};
      sel_out  <= cond ? hi : lo;
      hi_reg   <= hi;
      lo_reg   <= lo;
      cond_reg <= cond;
    end
  end

  assert property (@(posedge clk)
                    wide[7:4] == hi_reg && wide[3:0] == lo_reg
                    && ((cond_reg && sel_out == hi_reg)
                        || (!cond_reg && sel_out == lo_reg)));

endmodule
