// Tests a continuous assignment (`assign`).
// The wire `sum` is constrained to equal a + b.
// Property: sum != 30.  BMC must falsify with a=15, b=15.
module continuous_assign (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic [4:0] sum;
  assign sum = a + b;
  assert property (@(posedge clk) sum != 30);
endmodule
