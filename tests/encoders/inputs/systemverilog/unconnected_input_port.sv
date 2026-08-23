// A child instance's input port left explicitly unconnected (`.a()`).
// This encoder has no 4-state/high-Z representation, so an
// unconnected input is modeled as a free (undriven) input; `bout`,
// driven combinationally from `a`, should be freely choosable by BMC
// rather than stuck at a fixed value.
module child_unconn (input logic clk, input logic a, output logic b);
  assign b = a;
endmodule

module unconnected_input_port (input logic clk, input logic rst);

  logic bout;
  child_unconn c (.clk(clk), .a(), .b(bout));

  assert property (@(posedge clk) bout == 1'b0);

endmodule
