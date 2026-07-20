// Tests bit replication {{N{x}}}.
// Reset clears x; one cycle later, free input a = 2'b11 gives
// {{4{2'b11}}} = 8'b11111111 = 255, falsifying the property.
module replication (input logic clk, input logic rst, input logic [1:0] a);
  logic [7:0] x;
  always_ff @(posedge clk) begin
    if (rst) x <= 0;
    else x <= {4{a}};
  end
  assert property (@(posedge clk) x != 255);
endmodule
