// Out-of-scope: the `bind` directive, instantiating a checker-like
// module into another module's scope from outside its source.
module leaf2 (input logic clk, input logic rst);
  logic [3:0] count;
  always_ff @(posedge clk) begin
    if (rst) count <= '0;
    else count <= count + 1;
  end
endmodule

module checker_mod (input logic clk);
endmodule

module bind_directive (input logic clk, input logic rst);
  leaf2 l2 (
      clk,
      rst
  );
  bind leaf2 checker_mod cm (clk);
  assert property (@(posedge clk) 1'b1);
endmodule
