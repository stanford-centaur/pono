// goto_after_element.sv naming only the first of the two cycles it
// fires on. Must be refuted, or that fixture would hold for a
// counter that had simply stopped firing after the first match.
module goto_after_element_fails (
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
                   |-> ((b ##1 a[->2]) |-> (cnt == 5'd5)));

endmodule
