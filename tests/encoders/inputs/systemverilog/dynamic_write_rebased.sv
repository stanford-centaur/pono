// A runtime-indexed write into a packed range that does not start at
// zero. A declared index is a bit position only for an `[n:0]` range;
// for `[7:4]` it sits four above, so the write has to rebase before
// it can shift a mask into place. Without that, `r[4]` lands at bit
// position 4 -- off the top of a 4-bit value -- and the masked splice
// quietly drops it.
//
// `i` is driven to 4, so `r[i]` is the bottom bit of `r`.
module dynamic_write_rebased (
    input logic clk
);

  logic [7:4] r;
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

  assert property (@(posedge clk) started |-> (r[4] == 1'b1 && r == 4'b0001));

endmodule
