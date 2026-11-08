// Out-of-scope: the dynamic `string` type.
module string_type (input logic clk, input logic rst);
  string s;
  always_ff @(posedge clk) begin
    if (rst) s <= "";
    else s <= "x";
  end
  assert property (@(posedge clk) 1'b1);
endmodule
