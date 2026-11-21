// A counted repetition at `##0` after another element: the counting
// window opens on the prefix's own last cycle rather than the one
// after it. A record made there belongs to the count *before* that
// cycle's own occurrence, so it joins the register slot being
// shifted out rather than the fresh one -- and the window may be the
// one opening right now, which no register has seen yet.
//
// `b` is at cnt 3 and 4 and `a` at cnt 3, 5 and 8, so b@3 lands on
// an occurrence, which is exactly where `##0` and `##1` must part:
//
//   b ##0 a[->1] -> cnt 3 and 5      b ##1 a[->1] -> cnt 5
//   b ##0 a[->2] -> cnt 5 and 8      b ##1 a[->2] -> cnt 8
//
// Asserting all four together is what pins the difference; a
// `##0` quietly encoded as `##1` would fail the first and third.
module goto_after_element_overlap (
    input logic clk,
    input logic rst
);

  logic [4:0] cnt;
  logic a, b;

  always_ff @(posedge clk)
    cnt <= rst ? 5'd0 : (cnt == 5'd31 ? cnt : cnt + 5'd1);
  assign b = !rst && ((cnt == 5'd3) || (cnt == 5'd4));
  assign a = !rst && ((cnt == 5'd3) || (cnt == 5'd5) || (cnt == 5'd8));

  assert property (@(posedge clk) (!rst)
                   |-> (((b ##0 a[->1]) |-> (cnt == 5'd3 || cnt == 5'd5))
                        and ((b ##1 a[->1]) |-> (cnt == 5'd5))
                        and ((b ##0 a[->2]) |-> (cnt == 5'd5 || cnt == 5'd8))
                        and ((b ##1 a[->2]) |-> (cnt == 5'd8))));

endmodule
