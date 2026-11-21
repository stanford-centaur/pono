// The unknown half of xz_literal.sv: an x or z bit is unconstrained,
// not zero, so claiming a value for one must be refuted. Reading
// those bits as zeros -- the other plausible wrong answer, and the
// one that is silent -- would let this hold.
module xz_literal_fails (
    input logic clk,
    input logic rst
);

  logic [3:0] v;
  assign v = 4'b10xz;

  assert property (@(posedge clk) v[1] == 1'b0);

endmodule
