// The reading packed_array_ranges.sv rules out: the declared index
// taken as a bit offset counted up from zero, which puts `asc[0]` in
// the bottom byte and `off[2]`'s write two bytes above where it
// belongs. Exactly one of the two fixtures can hold, so refuting this
// one is what keeps the other from passing vacuously.
module packed_array_ranges_fails (
    input logic clk,
    input logic rst
);

  logic [0:3][7:0] asc;
  logic [5:2][7:0] off;

  always_ff @(posedge clk) begin
    asc <= 32'hAABBCCDD;

    off[5] <= 8'hAA;
    off[4] <= 8'hBB;
    off[3] <= 8'hCC;
    off[2] <= 8'hDD;
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (asc[0] == 8'hDD || off == 32'hDDCCBBAA));

endmodule
