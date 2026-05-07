// Tests bitwise operators (&, |, ^).
// Property: !(x == 15).  BMC must falsify with a=12, b=3
// (a&b=0, a|b=15, a^b=15 -> any of those reach 15).
module binary_bitwise (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic [3:0] x;
  always_ff @(posedge clk) x <= a | b;
  assert property (@(posedge clk) x != 15);
endmodule
