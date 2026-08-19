// A child instance's input port left explicitly unconnected (`.a()`).
// Per the LRM an unconnected input reads as high-impedance; since this
// encoder has no 4-state/high-Z representation, the closest correct
// model is to treat it as a completely free (undriven) input, the same
// as any other undriven signal elsewhere in this encoder -- so `bout`
// (driven combinationally from the child's now-free `a`) should be
// freely choosable by BMC, not stuck at some fixed value.
module child_unconn (input logic clk, input logic a, output logic b);
  assign b = a;
endmodule

module unconnected_input_port (input logic clk, input logic rst);

  logic bout;
  child_unconn c (.clk(clk), .a(), .b(bout));

  assert property (@(posedge clk) bout == 1'b0);

endmodule
