// `foreach` version of while_loop.sv, iterating directly over the
// bits of the packed vector `din` rather than a manually-managed
// index.
module foreach_loop (input logic clk, input logic rst, input logic [3:0] din);

  logic any_set, reg_any_set;

  always_comb begin
    any_set = 1'b0;
    foreach (din[i]) begin
      any_set |= din[i];
    end
  end

  always_ff @(posedge clk) begin
    if (rst) reg_any_set <= 1'b0;
    else reg_any_set <= any_set;
  end

  assert property (@(posedge clk) reg_any_set != 1'b1);

endmodule
