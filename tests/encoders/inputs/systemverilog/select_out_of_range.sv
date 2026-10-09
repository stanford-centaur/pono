// Reads that name bits the vector does not have. The LRM gives those
// bits X, so nothing may be proved about them -- but whatever the
// select does reach must still come through exactly.
//
// A constant index past the end used to reach the solver as an
// out-of-bounds Extract and abort; a dynamic one shifted in zeros
// instead, which made `p[i] == 0` provable for an i the array has no
// element for. The refutations live in select_out_of_range_fails.sv,
// since only the in-range halves can be asserted here.
module select_out_of_range (
    input logic clk,
    input logic [3:0][7:0] p,
    input logic [1:0] i,
    input logic [7:0] w
);

  logic [7:0] partial;
  assign partial = w[11:4];

  assert property (@(posedge clk) (p[2] == p[23:16] && p[i] == p[i]
                                   && partial[3:0] == w[7:4]));

endmodule
