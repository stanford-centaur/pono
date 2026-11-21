// Out-of-scope: constrained randomization (`rand`, `constraint`,
// `.randomize()`).
module randomize_constraint (input logic clk, input logic rst);
  class my_cls;
    rand int x;
    constraint c_range {
      x inside {[0 : 10]};
    }
  endclass
  my_cls obj;
  initial begin
    obj = new();
    void'(obj.randomize());
  end
  assert property (@(posedge clk) 1'b1);
endmodule
