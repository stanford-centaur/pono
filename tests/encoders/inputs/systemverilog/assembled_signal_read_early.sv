// reg_bus_slice.sv with the assertion moved above the instances that
// assemble `bus`. A module-level concurrent assertion is walked in
// the *first* pass -- it is a ProceduralBlockSymbol(Always) with no
// nonblocking targets and a non-edge-triggered body, so pass one's
// Always branch sends it to process_always_comb_once() -- which puts
// it in member order alongside the instances. Reading `bus` before
// they are walked throws "unknown symbol".
//
// Intended: the same verdict as reg_bus_slice.sv, which is this
// design with the assertion written last. Where a reader sits in the
// file should not change what the design means.
module child4_reg_early (input logic clk, input logic [3:0] in,
                         output logic [3:0] out);
  always_ff @(posedge clk) out <= in;
endmodule

module assembled_signal_read_early (input logic clk,
                                    input logic [3:0] a,
                                    input logic [3:0] b,
                                    input logic [3:0] c);

  logic [2:0][3:0] bus;
  logic started;

  initial started = 1'b0;

  assert property (@(posedge clk) started
                   |-> (bus[0] == $past(a) && bus[1] == $past(c)
                        && bus[2] == $past(b)));

  child4_reg_early u_lo  (.clk(clk), .in(a), .out(bus[0]));
  child4_reg_early u_hi  (.clk(clk), .in(b), .out(bus[2]));
  child4_reg_early u_mid (.clk(clk), .in(c), .out(bus[1]));

  always_ff @(posedge clk) started <= 1'b1;

endmodule
