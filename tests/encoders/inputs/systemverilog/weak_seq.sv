// `weak(seq)`: unlike strong_seq.sv, a `weak` sequence carries no
// obligation to ever match -- if `a` never fires, the sequence simply
// never starts, and the property holds vacuously forever.
module weak_seq (input logic clk, input logic rst, input logic a);

  assert property (@(posedge clk) weak(a ##1 1'b1));

endmodule
