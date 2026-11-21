// `assign`/`deassign`, the variable-only siblings of
// `force`/`release` (LRM 10.6.1). A procedural continuous assignment
// overrides procedural assignments to the variable until it is
// undone; applied at clock edges, its extent is in cycles.
module assign_deassign (
    input logic clk,
    input logic en,
    input logic [3:0] v
);

  logic [3:0] x;
  always_ff @(posedge clk) x <= x + 4'd1;

  always @(posedge clk) begin
    if (en) assign x = v;
    else deassign x;
  end

  assert property (@(posedge clk) en |=> x == $past(v));

endmodule
