// The same with a two-cycle delay, which moves where the window
// opens: from b@1 the window is [3,k] and still fills at k=5, but
// from b@4 it is [6,k], which never sees a second `a`. So the
// larger delay leaves exactly one firing cycle where `##1` had two
// -- a difference a counter keyed on the wrong cycle would miss.
//
// Also covers `[=n]`, which shares the registers but drops the
// "ends on an occurrence" half, so it keeps firing past cnt 5.
module goto_after_element_delay (
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
                   |-> (((b ##2 a[->2]) |-> (cnt == 5'd5))
                        and ((b ##1 a[=2]) |-> (cnt >= 5'd5))));

endmodule
