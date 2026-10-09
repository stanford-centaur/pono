// The other shape with no single internal symbol: an explicit port
// bound to a concatenation of several signals.
module explicit_port_concat (
    input logic clk,
    input logic a,
    input logic b,
    output .o({x, y})
);

  logic x, y;
  assign x = a;
  assign y = b;

  assert property (@(posedge clk) 1'b1);

endmodule
