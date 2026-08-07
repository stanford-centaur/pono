// TOP PRIORITY: `assume property` is supposed to constrain the trace
// so that `data` can never be anything but 0.  If the assumption is
// actually honored, `data != 7` holds forever (UNKNOWN).  `data` is
// otherwise a completely free input, so if the assumption is dropped,
// BMC can pick data == 7 immediately (falsifying at cycle 1) with
// nothing to stop it.
module assume_property (input logic clk, input logic rst, input logic [3:0] data);

  assume property (@(posedge clk) data == 4'd0);
  assert property (@(posedge clk) data != 4'd7);

endmodule
