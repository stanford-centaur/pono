// A combinational block that assigns a variable on some paths but
// not others. Real synthesis infers a latch there, and so does this
// now: the accumulated value already falls back to the symbol's own
// term where nothing wrote it, so committing it as a next-state
// update rather than a same-cycle equality is exactly that latch.
// Left as an equality the variable was simply free on the unwritten
// paths, which is sound but not what the hardware does.
//
// `y` is the latch; `w` is assigned unconditionally in the same
// block and must stay combinational, same cycle, with no added
// delay -- a mixed block is the case that decides this is a
// per-variable question rather than a per-block one.
module always_comb_latch (
    input logic clk,
    input logic rst,
    input logic en,
    input logic [7:0] d
);

  logic [7:0] w;
  logic [7:0] y;

  always_comb begin
    w = d + 8'd1;
    if (en) y = d;
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (w == d + 8'd1
                        && ($past(en) || y == $past(y))
                        && (!$past(en) || y == $past(d))));

endmodule
