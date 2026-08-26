// Out-of-scope: continuous `force`/`release` of a variable, which can
// override its always_ff-driven value from outside the block.
module force_release (input logic clk, input logic rst);
  logic [3:0] x;
  always_ff @(posedge clk) begin
    if (rst) x <= 4'd0;
    else x <= x + 4'd1;
  end
  initial begin
    force x = 4'd5;
    #1 release x;
  end
  assert property (@(posedge clk) 1'b1);
endmodule
