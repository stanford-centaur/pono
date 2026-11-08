// A child instance's output port connected to a concatenation of
// parent-side signals (`.sum({hi, lo})`), splitting the port's driven
// value across two separate parent nets: hi gets bit 3, lo gets bits
// 2:0, so {hi, lo} should be fixed at 4'b1010, not free.
module child_split (input logic clk, output logic [3:0] sum);
  assign sum = 4'b1010;
endmodule

module concat_output_port (input logic clk, input logic rst);

  logic hi;
  logic [2:0] lo;
  child_split c (.clk(clk), .sum({hi, lo}));

  assert property (@(posedge clk) {hi, lo} == 4'b1010);

endmodule
