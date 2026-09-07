// Companion to arith_compare.sv: same registers, but asserting the
// impossible direction of the inequality.  Already violated at
// cycle 1 by the reset values alone (0 < 0 is false).
module arith_compare_fails (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);

  logic [4:0] sum_reg;
  logic [3:0] a_prev;

  always_ff @(posedge clk) begin
    if (rst) begin
      sum_reg <= 5'd0;
      a_prev  <= 4'd0;
    end else begin
      sum_reg <= {1'b0, a} + {1'b0, b};
      a_prev  <= a;
    end
  end

  assert property (@(posedge clk) sum_reg < {1'b0, a_prev});

endmodule
