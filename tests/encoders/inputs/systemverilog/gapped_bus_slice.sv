// Three sibling instances each drive a different element of one shared
// packed array `bus`. u_lo writes element 0 first (a write starting
// at bit 0 but not covering the whole bus), then u_hi writes the
// non-adjacent element 2, then u_mid fills in element 1. Exercises
// process_continuous_assign_operand()'s full-width-write check, which
// must compare against the wire's declared width rather than mistake
// u_lo's partial first write for a full one.
module child4 (input logic [3:0] in, output logic [3:0] out);
  assign out = in;
endmodule

module gapped_bus_slice (input logic clk,
                          input logic [3:0] a,
                          input logic [3:0] b,
                          input logic [3:0] c);

  logic [2:0][3:0] bus;

  child4 u_lo  (.in(a), .out(bus[0]));
  child4 u_hi  (.in(b), .out(bus[2]));
  child4 u_mid (.in(c), .out(bus[1]));

  assert property (@(posedge clk) bus[0] == a && bus[1] == c && bus[2] == b);

endmodule
