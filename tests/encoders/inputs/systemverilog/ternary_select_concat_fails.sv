// Companion to ternary_select_concat.sv: asserts the concat/select
// inverse relationship does *not* hold.  Already violated at cycle 1
// (both sides are 0 there).
module ternary_select_concat_fails (
    input logic clk,
    input logic rst,
    input logic [3:0] hi,
    input logic [3:0] lo
);

  logic [7:0] wide;
  logic [3:0] hi_reg;

  always_ff @(posedge clk) begin
    if (rst) begin
      wide   <= 8'd0;
      hi_reg <= 4'd0;
    end else begin
      wide   <= {hi, lo};
      hi_reg <= hi;
    end
  end

  assert property (@(posedge clk) wide[7:4] != hi_reg);

endmodule
