// The dynamic half of the same story: a right shift feeds in zeros
// past the end of the vector, so this was provable. The LRM reads
// those bits as X, so it must be refuted.
module select_out_of_range_dynamic_fails (
    input logic clk,
    input logic [3:0][7:0] p,
    input logic [2:0] i
);

  assert property (@(posedge clk) (i >= 3'd4) |-> p[i] == 8'd0);

endmodule
