// Tests binary arithmetic operators (+, -, *).
// Reset clears the registers; one cycle later, free inputs a=3, b=2
// drive sum=5, diff=1, prod=6 simultaneously, falsifying the property.
module binary_arith (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic [3:0] sum, diff, prod;
  always_ff @(posedge clk) begin
    if (rst) begin
      sum  <= 0;
      diff <= 0;
      prod <= 0;
    end else begin
      sum  <= a + b;
      diff <= a - b;
      prod <= a * b;
    end
  end
  assert property (@(posedge clk) !(sum == 5 && diff == 1 && prod == 6));
endmodule
