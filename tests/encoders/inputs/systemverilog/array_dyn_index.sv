// Same packed array-of-vectors shape as array_const_index.sv, but both
// the write and the read use a *variable* element index (wsel/rsel are
// free inputs), exercising ElementSelectExpression with a
// non-constant index on both sides of an assignment.  Every cycle
// increments whichever element wsel points at; the array is all-zero
// at cycle 1 (the reset value), so the earliest an element can reach
// the max value 4'hF is after 15 *further* transitions all pointing
// wsel at the same element, i.e. at cycle 16 (one increment per
// transition, no wraparound below 16).
module array_dyn_index (
    input logic clk,
    input logic rst,
    input logic [1:0] wsel,
    input logic [1:0] rsel
);

  logic [3:0][3:0] mem;

  always_ff @(posedge clk) begin
    if (rst) mem <= '0;
    else mem[wsel] <= mem[wsel] + 4'd1;
  end

  assert property (@(posedge clk) mem[rsel] != 4'hF);

endmodule
