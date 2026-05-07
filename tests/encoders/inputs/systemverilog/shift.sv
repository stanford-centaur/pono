// Tests left-shift operator (<<).
// Reset clears x; one cycle later, free input a=3 makes x = 4'b1100 =
// 12 and falsifies the property.
module shift (input logic clk, input logic rst, input logic [3:0] a);
  logic [3:0] x;
  always_ff @(posedge clk) begin
    if (rst) x <= 0;
    else x <= a << 2;
  end
  assert property (@(posedge clk) x != 12);
endmodule
