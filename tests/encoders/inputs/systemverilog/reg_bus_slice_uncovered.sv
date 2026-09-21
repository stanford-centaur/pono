// A register bank assembled from sibling instances that leave a hole:
// `bus[1]` has no driver. The target is seeded with a free value and
// each contributor overwrites its own bits, so an uncovered bit keeps
// that free value and would drift every cycle with nothing saying so.
// Refused once every instance has been processed, which is the first
// point at which the hole is visible.
module child4_reg (
    input  logic       clk,
    input  logic [3:0] in,
    output logic [3:0] out
);

  always_ff @(posedge clk) out <= in;

endmodule

module reg_bus_slice_uncovered (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b
);

  logic [2:0][3:0] bus;

  child4_reg u_lo (.clk(clk), .in(a), .out(bus[0]));
  child4_reg u_hi (.clk(clk), .in(b), .out(bus[2]));

  assert property (@(posedge clk) 1'b1);

endmodule
