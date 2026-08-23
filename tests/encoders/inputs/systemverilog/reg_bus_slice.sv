// Three sibling instances, each an always_ff register, output-port
// aliased to a different element of one shared packed array `bus` --
// the register analogue of gapped_bus_slice.sv's wire-splicing case.
// Each instance's own piece only covers part of `bus`'s declared
// width, so declare_variables_internal() never creates a state var
// for `bus` at all: a mainstream RTL feature (a register bank driven
// element-by-element through an instance array) this encoder doesn't
// support, which must throw rather than silently produce a design
// with `bus` missing entirely.
module child4_reg (input logic clk, input logic [3:0] in, output logic [3:0] out);
  always_ff @(posedge clk) out <= in;
endmodule

module reg_bus_slice (input logic clk,
                      input logic [3:0] a,
                      input logic [3:0] b,
                      input logic [3:0] c);

  logic [2:0][3:0] bus;

  child4_reg u_lo  (.clk(clk), .in(a), .out(bus[0]));
  child4_reg u_hi  (.clk(clk), .in(b), .out(bus[2]));
  child4_reg u_mid (.clk(clk), .in(c), .out(bus[1]));

  // Deliberately doesn't reference `bus` here: without this fix, the
  // encoder would otherwise silently succeed with the three registers
  // simply absent from the model (confirmed empirically) rather than
  // failing at the point the gap is actually hit.
  assert property (@(posedge clk) 1'b1);

endmodule
