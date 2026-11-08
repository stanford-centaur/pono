// `$onehot`/`$onehot0`/`$isunknown`, checked as identities: exactly
// one bit set implies at most one bit set, for every `a`; and since
// Pono's SMT bitvector model has no X/Z state, `$isunknown` should
// always report false.
module onehot_isunknown (input logic clk, input logic rst, input logic [3:0] a);

  assert property (@(posedge clk) (!$onehot(a) || $onehot0(a)) && !$isunknown(a));

endmodule
