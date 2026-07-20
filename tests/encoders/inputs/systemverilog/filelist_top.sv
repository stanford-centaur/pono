// Top module for the --sv-filelist test: `child` lives in
// filelist_sub.sv, named by filelist.f rather than in this file.
//
// If the list-file's extra source is elaborated correctly, y == a + b
// at every cycle, so for a == 3, b == 5 we must have y == 8 != 0, and
// the assertion holds (BMC returns UNKNOWN within any finite bound).
//
// If the extra source were skipped, y would be left as a free input
// and BMC would falsify the property by picking y == 0.
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
