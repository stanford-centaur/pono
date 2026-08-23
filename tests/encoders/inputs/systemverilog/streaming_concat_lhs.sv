// A streaming concatenation used as an assignment target
// (`{>>{hi, lo}} <= a;`), which unpacks the RHS into hi/lo in stream
// order -- a distinct ExpressionKind (Streaming) from a plain
// concatenation-target LHS (ExpressionKind::Concatenation, already
// supported). resolve_lvalue() has no case for it, so this must throw
// a clear PonoException rather than silently dropping the write.
module streaming_concat_lhs (input logic clk, input logic [7:0] a);

  logic [3:0] hi, lo;

  always_ff @(posedge clk) begin
    { >> {hi, lo} } <= a;
  end

endmodule
