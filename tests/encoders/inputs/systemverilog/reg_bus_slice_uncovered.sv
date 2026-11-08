// A register bank assembled from sibling instances that leave a hole:
// `bus[1]` has no driver. An undriven signal is unconstrained, the
// same as an unconnected input port or an undriven net, so claiming a
// value for those bits must be refuted -- tying them to zero is a
// synthesis tool's choice, not the language's.
//
// Free every cycle rather than held, which is what a signal nothing
// drives does: a register would keep its value, and this is not one.
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

  assert property (@(posedge clk) bus[1] == 4'h0);

endmodule
