// Out-of-scope: the `real` (non-integral) type.
module real_type (input logic clk, input logic rst);
  real r;
  always_ff @(posedge clk) begin
    if (rst) r <= 0.0;
    else r <= r + 1.0;
  end
  assert property (@(posedge clk) 1'b1);
endmodule
