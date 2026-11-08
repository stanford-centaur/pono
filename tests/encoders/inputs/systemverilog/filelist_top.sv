// Top module for the --sv-filelist test: `child` lives in
// filelist_sub.sv, named by filelist.f rather than in this file.
// y == a + b at every cycle, so for a == 3, b == 5, y == 8 != 0
// always holds (BMC returns UNKNOWN). If the extra source were
// skipped, y would be a free input and BMC could pick y == 0,
// falsifying the property.
module top (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic [3:0] y;
  child c (.a(a), .b(b), .sum(y));
  assert property (@(posedge clk)
                   !(a == 4'd3 && b == 4'd5 && y == 4'd0));
endmodule
