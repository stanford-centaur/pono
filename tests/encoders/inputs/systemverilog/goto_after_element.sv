// A counted repetition *following* another element. The bare counter
// counts from the beginning of time, which is the wrong count here:
// the window opens where the preceding element ended. Writing A(t)
// for the occurrences through t, a window opening after cycle q
// holds exactly n of them at k when A(q) = A(k) - n -- so the
// question is whether the prefix ever ended at a cycle whose count
// was n less than now, answered by one register per count, shifted
// as the count rises. Only counts within n of the current one are
// ever asked about, so n + 1 registers suffice.
//
// Here `b` is at cnt 1 and 4, and `a` at cnt 3, 5 and 8, with the
// counter saturating so the pattern happens once. From b@1 the
// window [2,k] holds two a's at k=5; from b@4 the window [5,k] holds
// two at k=8. So `b ##1 a[->2]` fires at exactly cnt 5 and 8 --
// which the paired _fails fixtures pin by ruling out each alone.
module goto_after_element (
    input logic clk,
    input logic rst
);

  logic [4:0] cnt;
  logic a, b;

  always_ff @(posedge clk)
    cnt <= rst ? 5'd0 : (cnt == 5'd31 ? cnt : cnt + 5'd1);
  assign b = !rst && ((cnt == 5'd1) || (cnt == 5'd4));
  assign a = !rst && ((cnt == 5'd3) || (cnt == 5'd5) || (cnt == 5'd8));

  assert property (@(posedge clk) (!rst)
                   |-> ((b ##1 a[->2]) |-> (cnt == 5'd5 || cnt == 5'd8)));

endmodule
