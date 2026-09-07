// `restrict property` shares the exact same is_assumption/temporal
// handling as assume property in AssertionWalker::process_concurrent_
// assertion() -- see temporal_assume_property.sv.
module temporal_restrict_property (input logic clk, input logic a);
  restrict property (@(posedge clk) s_eventually a);
  assert property (@(posedge clk) 1'b1);
endmodule
