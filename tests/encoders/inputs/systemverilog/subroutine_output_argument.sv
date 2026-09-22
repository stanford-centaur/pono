// A function with an `output` formal, called from an expression. The
// write has to reach the caller, and an expression has no statement
// slot to deliver it in -- so the call has to be hoisted out, which
// is a statement-level change rather than an inliner one.
//
// Intended: provable, and the hand-hoisted equivalent is.
module subroutine_output_argument (input logic clk, input logic [7:0] d,
                                   output logic [7:0] q);

  function automatic [7:0] invert(input [7:0] x, output [7:0] dbl);
    dbl = x << 1;
    invert = x ^ 8'hFF;
  endfunction

  logic [7:0] t;
  always_comb q = invert(d, t);

  assert property (@(posedge clk) q == (d ^ 8'hFF) && t == (d << 1));

endmodule
