// De Morgan's law for bitwise operators (~(a&b) == ~a | ~b) and its
// logical-operator analogue ((a!=0)&&(b!=0) == !((a==0)||(b==0))),
// both registered from the same cycle's free a, b so the two sides of
// each identity are always fairly paired.  Genuinely holds for every
// a, b.  Paired with bitwise_logical_fails.sv.
module bitwise_logical (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);

  logic [3:0] demorgan_lhs, demorgan_rhs;
  logic       logic_lhs, logic_rhs;

  always_ff @(posedge clk) begin
    if (rst) begin
      demorgan_lhs <= 4'd0;
      demorgan_rhs <= 4'd0;
      logic_lhs    <= 1'b0;
      logic_rhs    <= 1'b0;
    end else begin
      demorgan_lhs <= ~(a & b);
      demorgan_rhs <= (~a) | (~b);
      logic_lhs    <= (a != 4'd0) && (b != 4'd0);
      logic_rhs    <= !((a == 4'd0) || (b == 4'd0));
    end
  end

  assert property (@(posedge clk)
                    demorgan_lhs == demorgan_rhs && logic_lhs == logic_rhs);

endmodule
