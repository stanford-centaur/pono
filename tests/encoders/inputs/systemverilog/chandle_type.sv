// Out-of-scope: the `chandle` (opaque C pointer) type.
module chandle_type (input logic clk, input logic rst);
  chandle h;
  always_ff @(posedge clk) begin
    if (rst) h <= null;
  end
  assert property (@(posedge clk) 1'b1);
endmodule
