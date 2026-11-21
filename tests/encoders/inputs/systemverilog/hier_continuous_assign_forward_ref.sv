// The same hierarchical drive as hier_continuous_assign.sv, with the
// assign written *above* the instance it reaches into. Declarations
// are walked in source order, so the target has no term yet when the
// assign is processed.
//
// Where a statement sits in the file does not change what the design
// means, so this holds exactly as its reordered twin does.
module hier_continuous_assign_fwd_sub;
  logic [1:0] x;
endmodule

module hier_continuous_assign_forward_ref (
    input logic clk
);

  assign u1.x[0] = 1'b1;
  hier_continuous_assign_fwd_sub u1 ();

  assert property (@(posedge clk) u1.x[0] == 1'b1);

endmodule
