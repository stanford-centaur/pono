// Arithmetic (+, -, *) and comparison (>=) composed as self-consistent
// invariants rather than a hand-picked falsifying trace: for any
// 4-bit a, b, (a+b) zero-extended to 5 bits is always >= a (b is
// unsigned, so it can only add); (a-b)+b == a mod 2^4 (subtraction is
// addition's inverse, exact through wraparound); and a*b == b*a
// (multiplication commutes).  All operands are registered from the
// *same* cycle's free inputs so each comparison is always fairly
// paired.  This holds for every a, b BMC could ever pick, so a
// genuine proof search (not just "no counterexample found yet")
// backs the UNKNOWN verdict.  Paired with arith_compare_fails.sv.
module arith_compare (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);

  logic [4:0] sum_reg;
  logic [3:0] a_prev, sub_plus_b, mul_ab, mul_ba;

  always_ff @(posedge clk) begin
    if (rst) begin
      sum_reg    <= 5'd0;
      a_prev     <= 4'd0;
      sub_plus_b <= 4'd0;
      mul_ab     <= 4'd0;
      mul_ba     <= 4'd0;
    end else begin
      sum_reg    <= {1'b0, a} + {1'b0, b};
      a_prev     <= a;
      sub_plus_b <= (a - b) + b;
      mul_ab     <= a * b;
      mul_ba     <= b * a;
    end
  end

  assert property (@(posedge clk)
                    sum_reg >= {1'b0, a_prev} && sub_plus_b == a_prev
                    && mul_ab == mul_ba);

endmodule
