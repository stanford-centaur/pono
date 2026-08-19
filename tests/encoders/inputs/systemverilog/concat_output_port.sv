// A child instance's output port connected to a concatenation of
// parent-side signals (`.sum({hi, lo})`), splitting the port's bits
// across two separate parent nets -- a legal, if less common, LRM
// connection shape. `find_lhs_base()`/`resolve_lvalue()` only handle a
// plain name/index/range/member as an output-port connection target;
// a concatenation should split the child's driven value across both
// halves (hi gets bit 3, lo gets bits 2:0), so hi/lo should end up
// fixed at 4'b1010's halves, not free.
module child_split (input logic clk, output logic [3:0] sum);
  assign sum = 4'b1010;
endmodule

module concat_output_port (input logic clk, input logic rst);

  logic hi;
  logic [2:0] lo;
  child_split c (.clk(clk), .sum({hi, lo}));

  assert property (@(posedge clk) {hi, lo} == 4'b1010);

endmodule
