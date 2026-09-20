// A temporary written only under a runtime condition has no value on
// the other path, and in a 2-valued model there is nothing to stand
// in for the X the LRM would give it.
module procedural_temp_cond_only (
    input logic clk,
    input logic c,
    input logic [7:0] a
);

  logic [7:0] y;

  always_ff @(posedge clk) begin
    logic [7:0] t;
    if (c) t = a;
    y <= t;
  end

  assert property (@(posedge clk) y == y);

endmodule
