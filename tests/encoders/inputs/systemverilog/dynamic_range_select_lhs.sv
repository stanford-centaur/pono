// A range-select lvalue with a non-constant base (`w[base +: 4]`,
// the indexed part-select form, which the LRM allows with a
// variable base as long as the width is constant). A fixed-width
// window at a runtime position is the same splice a runtime-indexed
// element select needs, only already counted in bits rather than
// elements -- so it takes the same dynamic-position path rather
// than a static bit range.
//
// The read side is the mirror image: shift the window down to the
// bottom and truncate.
//
// `started` gates the first cycle, where `$past` has no history and
// `w` is still free.
module dynamic_range_select_lhs (
    input logic clk,
    input logic [2:0] base,
    input logic [3:0] a
);

  logic [7:0] w;
  logic started;

  initial started = 1'b0;

  always_ff @(posedge clk) begin
    w[base+:4] <= a;
    started <= 1'b1;
  end

  // Only an in-range base: at 5 or more the window runs off the top
  // of `w`, which is a separate question from where it lands.
  assert property (@(posedge clk)
                   (started && $past(base) <= 3'd4)
                   |-> w[$past(base)+:4] == $past(a));

endmodule
