// The paired refutation for select_out_of_range.sv: element 7 of a
// four-element array does not exist, so it reads X and cannot equal
// any particular value. Before the bounds check this file could not
// even be encoded -- the Extract aborted inside the solver.
module select_out_of_range_fails (
    input logic clk,
    input logic [3:0][7:0] p
);

  assert property (@(posedge clk) p[7] == 8'd0);

endmodule
