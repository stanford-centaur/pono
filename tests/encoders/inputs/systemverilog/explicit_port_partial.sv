// An explicit port bound to a bit range rather than a whole signal.
// The port stands for four of `w`'s eight bits, and an internal
// symbol is the whole of one signal, so there is no term the two
// sides can share -- refused rather than bound to all of `w`.
module explicit_port_partial (
    input logic clk,
    input logic [7:0] d,
    output .o(w[3:0])
);

  logic [7:0] w;
  assign w = d;

  assert property (@(posedge clk) 1'b1);

endmodule
