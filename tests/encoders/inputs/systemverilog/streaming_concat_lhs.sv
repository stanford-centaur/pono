// A streaming concatenation used as an assignment target
// (`{>>{hi, lo}} <= a;`), which unpacks the RHS into hi/lo in stream
// order -- a distinct ExpressionKind (Streaming) from a plain
// concatenation-target LHS (ExpressionKind::Concatenation). It takes
// the same positional split, since `>>` re-orders nothing; what the
// stream form adds is that the source is consumed from its most
// significant end rather than truncated at the bottom.
//
// `started` gates the first cycle, where `$past` has no history.
module streaming_concat_lhs (
    input logic clk,
    input logic [7:0] a
);

  logic [3:0] hi, lo;
  logic started;

  initial started = 1'b0;

  always_ff @(posedge clk) begin
    {>>{hi, lo}} <= a;
    started <= 1'b1;
  end

  assert property (@(posedge clk) started
                   |-> hi == $past(a[7:4]) && lo == $past(a[3:0]));

endmodule
