// `force` and `release` applied at clock edges, so the override's
// extent is a number of cycles rather than an interval of simulation
// time. While `en` holds, `x` carries the forced value instead of
// the counter's.
module force_release_on_clock_edge (
    input logic clk,
    input logic en,
    input logic [3:0] v
);

  logic [3:0] x;
  always_ff @(posedge clk) x <= x + 4'd1;

  always @(posedge clk) begin
    if (en) force x = v;
    else release x;
  end

  assert property (@(posedge clk) en |=> x == $past(v));

endmodule
