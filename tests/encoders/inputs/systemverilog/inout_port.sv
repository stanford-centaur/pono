// An `inout` port. Only the child's drive is modelled, which is
// exactly right while nothing drives the net from outside -- the
// common single-direction-at-a-time use. What it cannot represent
// is the parent driving back, since resolving two drivers onto one
// net is not something this encoder has a notion of, so the
// approximation is logged rather than left silent.
//
// Here the child drives `io` high whenever `oe` is set, and that
// much must hold. While `oe` is low the child drives `z`, which
// this encoder models as an unconstrained bit, so nothing is
// claimed about `io` then -- inout_port_fails.sv pins that.
module inout_child (
    inout wire io,
    input logic oe
);
  assign io = oe ? 1'b1 : 1'bz;
endmodule

module inout_port (
    input logic clk,
    input logic rst,
    input logic oe
);

  wire io;
  inout_child u (.io(io), .oe(oe));

  assert property (@(posedge clk) oe |-> io);

endmodule
