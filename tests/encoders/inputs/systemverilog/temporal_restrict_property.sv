// `restrict property` shares the fairness-assumption path with
// `assume property` -- see temporal_assume_property.sv, of which
// this is the same design.
module temporal_restrict_property (input logic clk, input logic a);

  restrict property (@(posedge clk) s_eventually a);
  assert property (@(posedge clk) s_eventually a);

endmodule
