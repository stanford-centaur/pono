// Three sibling instances, each an always_ff register, output-port
// aliased to a different element of one shared packed array `bus` --
// the register analogue of gapped_bus_slice.sv's wire-splicing case.
//
// Each instance's own piece covers only part of `bus`, so aliasing
// the writes onto `bus` cannot work: all three would claim its
// next-state function and only one can have it. Each register keeps
// a state var of its own instead and splices its bits into `bus`,
// the way the wire case already assembles one. That the three of
// them between them cover `bus` is checked once every instance has
// been processed -- an uncovered bit would otherwise stay free.
//
// `started` gates the first cycle, where `$past` has no history.
module child4_reg (input logic clk, input logic [3:0] in, output logic [3:0] out);
  always_ff @(posedge clk) out <= in;
endmodule

module reg_bus_slice (input logic clk,
                      input logic [3:0] a,
                      input logic [3:0] b,
                      input logic [3:0] c);

  logic [2:0][3:0] bus;
  logic started;

  initial started = 1'b0;

  child4_reg u_lo  (.clk(clk), .in(a), .out(bus[0]));
  child4_reg u_hi  (.clk(clk), .in(b), .out(bus[2]));
  child4_reg u_mid (.clk(clk), .in(c), .out(bus[1]));

  always_ff @(posedge clk) started <= 1'b1;

  assert property (@(posedge clk) started
                   |-> (bus[0] == $past(a) && bus[1] == $past(c)
                        && bus[2] == $past(b)));

endmodule
