// Hierarchical design: `parent` instantiates `child`, which drives
// `parent.y` via its output port.
//
// If sub-modules are encoded correctly, y == a + b at every cycle,
// so for a == 3, b == 5 we must have y == 8 != 0, and the assertion
// holds (BMC returns UNKNOWN within any finite bound).
//
// If sub-modules are skipped (current encoder behavior), y is left
// as a free input and BMC will pick y == 0 with a == 3, b == 5,
// falsifying the property.
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
