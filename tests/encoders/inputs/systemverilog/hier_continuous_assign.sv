// A continuous assign driving a child instance's internal, non-port
// signal through a hierarchical dot-path, with the instance declared
// above the assign that reaches into it.
module hier_continuous_assign_sub;
  logic [1:0] x;
endmodule

module hier_continuous_assign (
    input logic clk
);

  hier_continuous_assign_sub u1 ();
  assign u1.x[0] = 1'b1;

  assert property (@(posedge clk) u1.x[0] == 1'b1);

endmodule
