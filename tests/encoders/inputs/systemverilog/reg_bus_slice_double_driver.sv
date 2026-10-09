// Two sibling instances whose registers both claim `bus[0]`. Each
// contributor to a spliced target writes its own bits, so two of
// them writing the same bits means two drivers on one signal, which
// nothing here resolves -- and silently letting one win would make
// the other instance's register vanish from the model.
module child4_reg (
    input  logic       clk,
    input  logic [3:0] in,
    output logic [3:0] out
);

  always_ff @(posedge clk) out <= in;

endmodule

module reg_bus_slice_double_driver (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b,
    input logic [3:0] c
);

  logic [2:0][3:0] bus;

  child4_reg u_lo (.clk(clk), .in(a), .out(bus[0]));
  child4_reg u_dup (.clk(clk), .in(b), .out(bus[0]));
  child4_reg u_mid (.clk(clk), .in(c), .out(bus[1]));
  child4_reg u_hi (.clk(clk), .in(b), .out(bus[2]));

  assert property (@(posedge clk) 1'b1);

endmodule
