// A latch is not a wire: `y` does not take `d` in the same cycle
// the enable is high, it takes it on the next edge. Must be
// refuted, or the latch would be indistinguishable from the
// combinational reading it replaces.
module always_comb_latch_fails (
    input logic clk,
    input logic rst,
    input logic en,
    input logic [7:0] d
);

  logic [7:0] y;

  always_comb begin
    if (en) y = d;
  end

  assert property (@(posedge clk) en |-> y == d);

endmodule
