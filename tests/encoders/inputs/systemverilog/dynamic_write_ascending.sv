// The ascending twin of dynamic_write_rebased: for `[4:7]` the
// leftmost declared index is the most significant bit, so bit
// position and declared index count in opposite directions and the
// rebasing is a subtraction from the top rather than from the
// bottom. Getting the direction wrong puts `r[4]` at the bottom of
// `r` instead of the top, which `r == 4'b1000` catches.
module dynamic_write_ascending (
    input logic clk
);

  logic [4:7] r;
  logic [2:0] i;
  logic started;

  initial begin
    started = 1'b0;
    r = 4'b0000;
    i = 3'd4;
  end

  always_ff @(posedge clk) begin
    i <= 3'd4;
    r[i] <= 1'b1;
    started <= 1'b1;
  end

  assert property (@(posedge clk) started |-> (r[4] == 1'b1 && r == 4'b1000));

endmodule
