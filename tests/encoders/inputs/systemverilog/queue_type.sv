// Out-of-scope: a queue (`int q[$];`).
module queue_type (input logic clk, input logic rst);
  int q[$];
  always_ff @(posedge clk) begin
    if (rst) q.push_back(1);
  end
  assert property (@(posedge clk) 1'b1);
endmodule
