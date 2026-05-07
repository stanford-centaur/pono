// Tests if/else inside an always_ff block.
// When sel = 1, x is set to 5; when sel = 0, x is set to 3.
// Property: x != 5.  BMC must falsify within 2 steps (sel=1).
module if_else (input logic clk, input logic sel);
  logic [3:0] x;
  always_ff @(posedge clk) begin
    if (sel) x <= 5;
    else x <= 3;
  end
  assert property (@(posedge clk) x != 5);
endmodule
