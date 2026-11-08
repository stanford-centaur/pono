// Four unrelated-looking identities bundled into one register-path
// invariant: a shift-left-by-1 equals multiply-by-2 mod 2^4 (exact,
// even through overflow), double bitwise-negation is the identity,
// unary minus is an additive inverse mod 2^4 (also exact through
// overflow), and reduction-OR agrees with an explicit nonzero check.
// All derived from the same cycle's free `a`, so genuinely holds for
// every a.  Paired with shift_unary_reduction_fails.sv.
module shift_unary_reduction (input logic clk, input logic rst, input logic [3:0] a);

  logic [3:0] shift_lhs, shift_rhs, dbl_not, neg_a, a_reg;
  logic       red_or, or_check;

  always_ff @(posedge clk) begin
    if (rst) begin
      shift_lhs <= 4'd0;
      shift_rhs <= 4'd0;
      dbl_not   <= 4'd0;
      neg_a     <= 4'd0;
      a_reg     <= 4'd0;
      red_or    <= 1'b0;
      or_check  <= 1'b0;
    end else begin
      shift_lhs <= a << 1;
      shift_rhs <= a * 2;
      dbl_not   <= ~(~a);
      neg_a     <= -a;
      a_reg     <= a;
      red_or    <= |a;
      or_check  <= (a != 4'd0);
    end
  end

  assert property (@(posedge clk)
                    shift_lhs == shift_rhs && dbl_not == a_reg
                    && (neg_a + a_reg) == 4'd0 && red_or == or_check);

endmodule
