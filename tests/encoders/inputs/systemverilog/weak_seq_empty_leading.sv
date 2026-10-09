// A sequence whose *leading* repetition can match emptily. The
// sequence itself never does -- the `##1 y` is always there -- but
// the leading `x` need not hold for an attempt to be under way, so
// there is no single cycle that marks where one began. weak() used
// to need exactly that and refused the shape; it now checks the
// attempt that began S cycles ago, and every cycle begins one.
//
// Three alternatives, of spans 0, 1 and 2, so S is 2 and the
// earliest cycle at which the attempt starting at cycle 0 is
// decided -- hence the earliest possible refutation -- is cycle 2.
module weak_seq_empty_leading (
    input logic clk,
    input logic x,
    input logic y
);

  assert property (@(posedge clk) weak(x[*0:2] ##1 y));

endmodule
