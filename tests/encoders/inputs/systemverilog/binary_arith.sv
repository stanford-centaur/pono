// Tests binary arithmetic operators (+, -, *).
// Property: !(sum == 5 && diff == 1 && prod == 6).
// BMC must falsify with a=3, b=2 (sum=5, diff=1, prod=6).
module binary_arith (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic [3:0] sum, diff, prod;
  always_ff @(posedge clk) begin
    sum  <= a + b;
    diff <= a - b;
    prod <= a * b;
  end
  assert property (@(posedge clk) !(sum == 5 && diff == 1 && prod == 6));
endmodule
