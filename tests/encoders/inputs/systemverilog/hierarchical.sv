// Hierarchical design: `parent` instantiates `child`, which drives
// `parent.y` via its output port. y == a + b at every cycle, so for
// a == 3, b == 5, y == 8 != 0 always holds (BMC returns UNKNOWN
// within any finite bound). If sub-modules were left unencoded, y
// would be a free input and BMC would falsify the property by
// picking y == 0.
module child (
    input  logic [3:0] a,
    input  logic [3:0] b,
    output logic [3:0] sum
);
  assign sum = a + b;
endmodule

module parent (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b
);
  logic [3:0] y;
  child c (.a(a), .b(b), .sum(y));
  assert property (@(posedge clk)
                   !(a == 4'd3 && b == 4'd5 && y == 4'd0));
endmodule
