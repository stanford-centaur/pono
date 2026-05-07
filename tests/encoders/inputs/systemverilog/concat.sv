// Tests bit concatenation {a, b}.
// Property: x != 171.  BMC must falsify with a=10, b=11
// (since {4'd10, 4'd11} = 8'b10101011 = 171).
module concat (input logic clk, input logic [3:0] a, input logic [3:0] b);
  logic [7:0] x;
  always_ff @(posedge clk) x <= {a, b};
  assert property (@(posedge clk) x != 171);
endmodule
