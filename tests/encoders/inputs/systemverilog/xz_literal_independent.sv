// Each unknown bit is its own free value, so two `x` literals are not
// forced to agree -- the loosest reading, and so the soundest. Must
// be refuted.
module xz_literal_independent (
    input logic clk,
    input logic rst
);

  logic a, b;
  assign a = 1'bx;
  assign b = 1'bx;

  assert property (@(posedge clk) a == b);

endmodule
