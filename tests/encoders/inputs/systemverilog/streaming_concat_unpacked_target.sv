// A streaming concatenation whose target is an unpacked array, which
// is walked in `foreach` order (LRM 11.4.14.1) rather than as a flat
// bit range. Fixed-size unpacked arrays are modelled here, so the
// elements the stream lands in are known.
//
// `>>` re-orders nothing, so the source's top bits go to the first
// element.
module streaming_concat_unpacked_target (
    input logic clk,
    input logic [7:0] a
);

  logic [3:0] mem[0:1];
  always_ff @(posedge clk) {>>{mem}} <= a;

  logic started;
  initial started = 1'b0;
  always_ff @(posedge clk) started <= 1'b1;

  assert property (@(posedge clk) started
                   |-> (mem[0] == $past(a[7:4]) && mem[1] == $past(a[3:0])));

endmodule
