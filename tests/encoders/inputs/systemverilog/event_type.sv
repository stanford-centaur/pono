// Out-of-scope: the `event` type and the `->` trigger statement.
module event_type (input logic clk, input logic rst);
  event e;
  always_ff @(posedge clk) begin
    if (rst) -> e;
  end
  assert property (@(posedge clk) 1'b1);
endmodule
