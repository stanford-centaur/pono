// Hierarchical value reference: the parent reads the child's
// internal counter via a dotted name (c.cnt) instead of through a
// port.  The assertion `c.cnt != 5` falsifies once the counter
// reaches 5, which (with a one-cycle reset) happens at cycle 6.
module child_counter (input logic clk, input logic rst);
  logic [3:0] cnt;
  always @(posedge clk) cnt <= rst ? '0 : cnt + 1;
endmodule

module hierarchical_value_top (input logic clk, input logic rst);
  child_counter c (.clk(clk), .rst(rst));
  assert property (@(posedge clk) c.cnt != 4'd5);
endmodule
