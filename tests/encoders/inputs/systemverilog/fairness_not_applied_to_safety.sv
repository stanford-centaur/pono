// Where a fairness assumption stops. `dead` latches and from then
// on holds `a` low, so no trace through it satisfies `s_eventually a`
// infinitely often -- under an infinite-trace reading `dead` is
// unreachable, and fairness_applied_to_liveness.sv proves exactly
// that of this same design.
//
// A safety counterexample is a finite trace, though, and no finite
// trace contradicts a fairness constraint: LRM 1800-2009 F.5.3.1
// makes an assert hold on every *feasible* word, and F.5.3.2 makes
// `s_eventually a` weakly satisfied by every finite word, so every
// finite prefix is feasible. This is also the split BTOR2 makes,
// where `fair` joins a justice property's conditions and a `bad`
// state is checked without it. So this is refuted, at the one step
// it takes `go` to set `dead`.
module fairness_not_applied_to_safety (
    input logic clk,
    input logic a,
    input logic go
);

  logic dead;
  initial dead = 1'b0;
  always_ff @(posedge clk) if (go) dead <= 1'b1;

  // Once dead, the fairness condition can never be met again.
  assume property (@(posedge clk) dead |-> !a);
  assume property (@(posedge clk) s_eventually a);

  assert property (@(posedge clk) !dead);

endmodule
