// reg_bus_slice.sv with the assertion moved above the instances that
// assemble `bus`, so the reader is walked before its contributors.
// Module-level concurrent assertions are walked in the same pass as
// instances and in member order, which is what makes the two
// orderings different designs to the encoder rather than the same
// one. They describe the same hardware, so they hold the same way.
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
