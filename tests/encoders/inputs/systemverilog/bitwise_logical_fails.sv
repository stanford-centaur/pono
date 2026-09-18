// Companion to bitwise_logical.sv: asserts the two De Morgan sides
// *differ*.  Already violated at cycle 1 (both sides are 0 there).
module bitwise_logical_fails (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);

  logic [3:0] demorgan_lhs, demorgan_rhs;

  always_ff @(posedge clk) begin
    if (rst) begin
      demorgan_lhs <= 4'd0;
      demorgan_rhs <= 4'd0;
    end else begin
      demorgan_lhs <= ~(a & b);
      demorgan_rhs <= (~a) | (~b);
    end
  end

  assert property (@(posedge clk) demorgan_lhs != demorgan_rhs);

endmodule
