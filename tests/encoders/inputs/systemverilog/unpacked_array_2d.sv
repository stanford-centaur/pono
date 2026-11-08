// Unpacked arrays with more than one dimension. The sort is an array
// of arrays, so each index is a Select on the way down and a Store on
// the way back out; only the innermost write replaces a value, and
// every enclosing dimension is rebuilt around it.
//
// The dimensions are deliberately unequal, and the assertions name
// neighbours in both of them, so an index applied to the wrong
// dimension either disturbs a neighbour or lands outside the array
// rather than quietly still working.
//
// `grid` is written at a runtime index; `fixed`'s neighbours pin the
// isolation; `deep` has three dimensions; `sub` writes a bit range
// inside a two-dimensional element; `row` is replaced whole, once
// from an assignment pattern and once from another array.
module unpacked_array_2d (
    input logic clk,
    input logic rst,
    input logic [1:0] i,
    input logic [2:0] j,
    input logic [7:0] d,
    input logic [3:0] v
);

  logic [7:0] grid[0:2][0:7];
  logic [7:0] fixed[0:3][0:3];
  logic [7:0] deep[0:1][0:1][0:3];
  logic [7:0] sub[0:3][0:3];
  logic [7:0] row[0:3][0:3];
  logic [7:0] src[0:3];

  always_ff @(posedge clk) begin
    if (rst) begin
      fixed <= '{default: 0};
      sub   <= '{default: 0};
    end else begin
      fixed[1][1] <= d;
      sub[2][3][3:0] <= v;
    end
    // Unconditional, so these hold one cycle after reset releases
    // without needing the reset value carried through.
    grid[2][j] <= d;
    deep[1][0][i] <= d;
    src <= '{default: 3};
    row[1] <= '{default: 7};
    row[2] <= src;
  end

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> (grid[2][$past(j)] == $past(d)
                        && fixed[1][1] == $past(d)
                        && fixed[1][2] == $past(fixed[1][2])
                        && fixed[2][1] == $past(fixed[2][1])
                        && deep[1][0][$past(i)] == $past(d)
                        && sub[2][3][3:0] == $past(v)
                        && row[1][i] == 8'd7
                        && row[2][i] == 8'd3));

endmodule
