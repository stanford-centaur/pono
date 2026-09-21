// The chain weak_seq_and.sv's `and` has to reduce to. Kept as its
// own fixture so the pair is the assertion rather than a depth
// written down in the test.
module weak_seq_and_chain (
    input logic clk,
    input logic x,
    input logic y,
    input logic z
);

  assert property (@(posedge clk) weak(x ##1 y ##1 z));

endmodule
