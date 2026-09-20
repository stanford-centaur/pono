// A repetition with a zero lower bound also matches emptily, and the
// LRM absorbs one cycle of the following delay when it does:
// `(empty ##1 b)` is `##0 b`, i.e. `b` on its own. So
// `b ##1 a[*0:2]` completes at cycle k when
//
//   b(k)                          -- the empty match, (b ##1 empty) = b
//   or (b(k-1) && a(k))           -- a[*1]
//   or (b(k-2) && a(k-1) && a(k)) -- a[*2]
//
// and `ok` is that expansion written out by hand. A trailing
// repetition is what makes this discriminating: a *leading* empty
// match drops every mention of `a`, so the expansion would collapse
// to `b` and could not tell a correct encoding from one that matches
// too often.
//
// Paired with empty_repeat_match_fails.sv, which shifts `ok` by a
// cycle so that this fixture cannot hold vacuously.
module empty_repeat_match (
    input logic clk,
    input logic rst,
    input logic a,
    input logic b
);

  logic ok;
  assign ok = b || ($past(b) && a) || ($past(b, 2) && $past(a) && a);

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> ((b ##1 a[*0:2]) |-> ok));

endmodule
