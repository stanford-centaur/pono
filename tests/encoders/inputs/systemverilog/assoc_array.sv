// Out-of-scope: an associative array (`int m[string];`).
module assoc_array (input logic clk, input logic rst);
  int m[string];
  always_ff @(posedge clk) begin
    if (rst) m["k"] = 1;
  end
  assert property (@(posedge clk) 1'b1);
endmodule
