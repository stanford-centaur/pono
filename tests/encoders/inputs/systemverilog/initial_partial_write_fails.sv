// The paired refutation for initial_partial_write.sv: x is provably
// 8'h08, so this must be refuted. Before the initial block's writes
// composed, they contradicted each other and nothing over the design
// could be refuted at all.
module initial_partial_write_fails (
    input logic clk
);

  logic [7:0] x;

  initial begin
    x = 8'd0;
    x[3] = 1'b1;
  end

  always_ff @(posedge clk) x <= x;

  assert property (@(posedge clk) x == 8'hFF);

endmodule
