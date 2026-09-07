// Procedural immediate `cover` (distinct from concurrent `cover
// property (...)`, covered by cover_property.sv). Same reachability-
// duality treatment (see sva.cpp): the negated condition is pushed
// into propvec() like a plain assertion, so a counterexample to "never
// covered" is exactly "the cover point was reached". `data` is free,
// so the cover point (data == 5) is reachable at the earliest possible
// cycle.
module immediate_cover (input logic clk, input logic rst, input logic [3:0] data);

  always_comb begin
    cover (data == 4'd5);
  end

endmodule
