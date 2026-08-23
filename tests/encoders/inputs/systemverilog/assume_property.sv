// `assume property` constrains the trace so `data` can never be
// anything but 0.  With the assumption honored, `data != 7` holds
// forever (UNKNOWN).  `data` is otherwise a completely free input, so
// if the assumption were dropped, BMC could pick data == 7 immediately
// (falsifying at cycle 1) with nothing to stop it.
module assume_property (input logic clk, input logic rst, input logic [3:0] data);

  assume property (@(posedge clk) data == 4'd0);
  assert property (@(posedge clk) data != 4'd7);

endmodule
