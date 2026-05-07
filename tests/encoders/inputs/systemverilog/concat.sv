// Tests bit concatenation {a, b}.
// Reset clears x; one cycle later, free inputs a=10, b=11 give
// {4'd10, 4'd11} = 8'b10101011 = 171, falsifying the property.
module concat (
    input logic clk,
    input logic rst,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic [7:0] x;
  always_ff @(posedge clk) begin
    if (rst) x <= 0;
    else x <= {a, b};
  end
  assert property (@(posedge clk) x != 171);
endmodule
