// Companion to shift_unary_reduction.sv: asserts the shift/multiply
// identity does *not* hold.  Already violated at cycle 1 (both 0).
module shift_unary_reduction_fails (input logic clk, input logic rst, input logic [3:0] a);

  logic [3:0] shift_lhs, shift_rhs;

  always_ff @(posedge clk) begin
    if (rst) begin
      shift_lhs <= 4'd0;
      shift_rhs <= 4'd0;
    end else begin
      shift_lhs <= a << 1;
      shift_rhs <= a * 2;
    end
  end

  assert property (@(posedge clk) shift_lhs != shift_rhs);

endmodule
