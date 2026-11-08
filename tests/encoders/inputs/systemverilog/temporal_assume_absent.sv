// temporal_assume_property.sv without the assumption. `a` is free
// and can stay low forever, so the eventuality is refuted.
module temporal_assume_absent (input logic clk, input logic a);

  assert property (@(posedge clk) s_eventually a);

endmodule
