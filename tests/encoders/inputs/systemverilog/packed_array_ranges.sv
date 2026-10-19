// Packed arrays whose declared range is not the conventional `[n:0]`.
// The left bound is the most significant end whichever way the range
// runs, so element 0 of `[0:3]` is the *top* byte, and `[5:2]`'s
// lowest element is 2, not 0. Using the declared index directly as a
// bit offset -- right only for `[n:0]` -- reads and writes a real but
// wrong element, so this fixture pins each shape against a byte
// pattern whose bytes are all distinct.
//
// `asc` is written whole and read element-wise, `off` and `neg` the
// other way round, so the read and write paths are each covered by at
// least one non-conventional range. `sel` is free, making the
// dynamic-index read a case split over every element rather than a
// single one.
module packed_array_ranges (
    input logic clk,
    input logic rst,
    input logic [1:0] sel
);

  logic [0:3][7:0] asc;
  logic [5:2][7:0] off;
  logic [1:-2][7:0] neg;
  logic [7:0] dyn;

  always_ff @(posedge clk) begin
    asc <= 32'hAABBCCDD;
    dyn <= asc[sel];

    off[5] <= 8'hAA;
    off[4] <= 8'hBB;
    off[3] <= 8'hCC;
    off[2] <= 8'hDD;

    neg[1]  <= 8'hAA;
    neg[0]  <= 8'hBB;
    neg[-1] <= 8'hCC;
    neg[-2] <= 8'hDD;
  end

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> (asc[0] == 8'hAA && asc[1] == 8'hBB
                        && asc[2] == 8'hCC && asc[3] == 8'hDD
                        && off == 32'hAABBCCDD
                        && neg == 32'hAABBCCDD
                        && dyn == ($past(sel) == 2'd0 ? 8'hAA
                                   : $past(sel) == 2'd1 ? 8'hBB
                                   : $past(sel) == 2'd2 ? 8'hCC : 8'hDD)));

endmodule
