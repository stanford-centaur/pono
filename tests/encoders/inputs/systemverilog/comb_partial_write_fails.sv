// The negative half of comb_partial_write.sv: p is provably 8'h08, so
// this must be refuted. If the block's writes became separate
// conflicting constraints, the design would be unsatisfiable and this
// would go unrefuted instead.
module comb_partial_write_fails (
    input logic clk
);

  logic [7:0] p;

  always_comb begin
    p = 8'd0;
    p[3] = 1'b1;
  end

  assert property (@(posedge clk) p == 8'hFF);

endmodule
