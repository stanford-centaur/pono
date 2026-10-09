// A plain `always` block is a register block or combinational logic
// depending only on its event control, and nothing else distinguishes
// the two when every write uses blocking `=`. Both spellings appear
// here: the edge-sensitive one must become a flop (x lags v by a
// cycle), the level-sensitive one must stay combinational (y tracks v
// in the same cycle).
module blocking_always_edge (
    input logic clk,
    input logic rst,
    input logic [7:0] v
);

  logic [7:0] x;
  logic [7:0] y;

  always @(posedge clk) x = v;
  always @(*) y = v + 8'd1;

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (x == $past(v) && y == v + 8'd1));

endmodule
