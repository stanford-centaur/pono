// The `##0` window read as if it opened a cycle later: cnt 3 is a
// match of `b ##0 a[->1]` and must not be excluded. Refuting this
// is what shows the overlap really is handled rather than folded
// into the `##1` case.
module goto_after_element_overlap_fails (
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
                   |-> ((b ##0 a[->1]) |-> (cnt == 5'd5)));

endmodule
