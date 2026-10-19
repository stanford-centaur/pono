// weak() of `or`: the match set is the union of the two operands'
// (LRM 16.9.7), so an attempt completes if either alternative does,
// and the check is anchored at the longer one's span -- 2 here, so
// the attempt that began 2 cycles ago is satisfied by `y` at either
// of the two cycles since.
module weak_seq_or (
    input logic clk,
    input logic x,
    input logic y
);

  assert property (@(posedge clk) weak((x ##1 y) or (x ##2 y)));

endmodule
