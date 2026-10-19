// The other half of inout_port.sv: with the child not driving, `io`
// is a `z` and so unconstrained, not low. Claiming it is low must be
// refuted -- reading an undriven net as zero is the silent wrong
// answer this rules out.
module inout_child_f (
    inout wire io,
    input logic oe
);
  assign io = oe ? 1'b1 : 1'bz;
endmodule

module inout_port_fails (
    input logic clk,
    input logic rst,
    input logic oe
);

  wire io;
  inout_child_f u (.io(io), .oe(oe));

  assert property (@(posedge clk) !oe |-> !io);

endmodule
