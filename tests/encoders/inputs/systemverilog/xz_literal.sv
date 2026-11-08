// An `x` or `z` bit in a literal. This encoder models X as an
// unconstrained value throughout (`===` collapses to `==`,
// `$isunknown` is constantly false, an unreset register is free at
// time 0), so an unknown bit in a literal becomes a fresh
// unconstrained bit and the known bits around it are pinned as
// written. The alternative was handing the solver a decimal string
// with an `x` in it, which aborts it rather than raising a
// PonoException.
//
// Here `v[3]` and `v[2]` are given, so they must hold; `v[1:0]` are
// not, and xz_literal_fails.sv pins that they are genuinely free
// rather than quietly zero.
module xz_literal (
    input logic clk,
    input logic rst
);

  logic [3:0] v;
  assign v = 4'b10xz;

  assert property (@(posedge clk) v[3] == 1'b1 && v[2] == 1'b0);

endmodule
