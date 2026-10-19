// Out-of-scope: a class handle variable, declared and constructed.
module class_decl (input logic clk, input logic rst);
  class my_cls;
    int x;
  endclass
  my_cls c;
  always_ff @(posedge clk) begin
    if (rst) c = new();
  end
  assert property (@(posedge clk) 1'b1);
endmodule
