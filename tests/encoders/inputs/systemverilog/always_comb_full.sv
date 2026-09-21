// The control for always_comb_latch.sv: assigned on every path, so
// no latch is inferred and the value is available in the same
// cycle. Treating every combinational target as a latch would add a
// cycle of delay here and this would fail.
module always_comb_full (
    input logic clk,
    input logic rst,
    input logic en,
    input logic [7:0] d
);

  logic [7:0] y;

  always_comb begin
    if (en) y = d;
    else y = 8'd0;
  end

  assert property (@(posedge clk) y == (en ? d : 8'd0));

endmodule
