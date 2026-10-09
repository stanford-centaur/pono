// Streaming an unpacked array into or out of a stream requires
// walking its elements in `foreach` order (LRM 11.4.14.1), which this
// encoder's positional bit-splicing has no representation for. The
// target must be integral, and saying so is what keeps the write from
// being read as a zero-width one.
module streaming_concat_unpacked_target (
    input logic clk,
    input logic [7:0] a
);

  logic [3:0] mem[0:1];
  logic [3:0] hi;

  always_ff @(posedge clk) {>>{hi, mem}} <= {a, a};

  assert property (@(posedge clk) 1'b1);

endmodule
