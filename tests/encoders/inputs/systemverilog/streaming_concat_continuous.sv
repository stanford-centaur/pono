// A streaming concatenation as a continuous-assignment target, and
// the wider-source rule that goes with it: a stream is consumed from
// its most significant end, so a 12-bit source fills an 8-bit target
// from a[11:4] and leaves the bottom four bits unread -- the opposite
// of the truncation every other assignment target gets.
module streaming_concat_continuous (
    input logic clk,
    input logic [11:0] a
);

  logic [3:0] hi, lo;

  assign {>>{hi, lo}} = a;

  assert property (@(posedge clk) hi == a[11:8] && lo == a[7:4]);

endmodule
