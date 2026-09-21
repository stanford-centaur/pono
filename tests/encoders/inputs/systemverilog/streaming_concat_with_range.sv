// A `with` range inside a streaming concatenation selects how much of
// an unpacked array takes part (LRM 11.4.14.4). It exists to size
// dynamically sized data, which is not part of this encoder's
// subject matter, so it is refused rather than ignored.
module streaming_concat_with_range (
    input logic clk,
    input logic [7:0] a
);

  logic [3:0] mem[0:1];

  always_ff @(posedge clk) {>>{mem with [0:1]}} <= a;

  assert property (@(posedge clk) 1'b1);

endmodule
