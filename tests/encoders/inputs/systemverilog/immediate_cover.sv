// Procedural immediate `cover` (distinct from concurrent `cover
// property (...)`, covered by cover_property.sv). Mirrors that
// fixture's assumed contract: a correct implementation should expose a
// cover goal's *reachability* through the same propvec()-of-safety-
// properties interface used for assert/assume, by pushing the negated
// condition (finding a counterexample to "never covered" is exactly
// "the cover point was reached"). `data` is free, so the cover point
// (data == 5) is reachable at the earliest possible cycle.
module immediate_cover (input logic clk, input logic rst, input logic [3:0] data);

  always_comb begin
    cover (data == 4'd5);
  end

endmodule
