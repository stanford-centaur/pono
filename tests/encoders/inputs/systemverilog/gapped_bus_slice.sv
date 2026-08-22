// Three sibling instances each drive a different element of one shared
// packed array `bus`, declared so the *first*-processed instance
// (u_lo) writes element 0 -- a write starting at bit 0 -- while the
// *second*-processed instance (u_hi) writes the non-adjacent element 2,
// with element 1 filled in last by u_mid. This exercises
// process_continuous_assign()'s "is this the wire's full-width write"
// check on a wire's very first write, for a write that starts at bit 0
// but does not cover the whole bus.
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
