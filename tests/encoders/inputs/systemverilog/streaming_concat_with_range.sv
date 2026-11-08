// A `with` range inside a streaming concatenation (LRM 11.4.14.4),
// selecting how much of an unpacked array takes part. The range here
// is constant and the array fixed-size, so which elements it names is
// settled at elaboration.
//
// `mem with [0:1]` covers both elements, so this is the whole array
// in `foreach` order.
module streaming_concat_with_range (
    input logic clk,
    input logic [7:0] a
);

  logic [3:0] mem[0:1];
  always_ff @(posedge clk) {>>{mem with [0:1]}} <= a;

  logic started;
  initial started = 1'b0;
  always_ff @(posedge clk) started <= 1'b1;

  assert property (@(posedge clk) started
                   |-> (mem[0] == $past(a[7:4]) && mem[1] == $past(a[3:0])));

endmodule
