// Tests if/else inside an always_ff block.
// Reset clears x; one cycle later, sel=1 sets x to 5 and falsifies
// the property.
module if_else (input logic clk, input logic rst, input logic sel);
  logic [3:0] x;
  always_ff @(posedge clk) begin
    if (rst) x <= 0;
    else if (sel) x <= 5;
    else x <= 3;
  end
  assert property (@(posedge clk) x != 5);
endmodule
