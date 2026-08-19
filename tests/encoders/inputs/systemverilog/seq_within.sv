// Sequence `within`: per the LRM, `s1 within s2` matches over the same
// span as a match of s2, provided s1 matches ending somewhere inside
// that span. For `b ##1 c` (spanning [i, i+1]) and single-cycle `a`,
// that reduces to: b@i AND c@i+1 AND (a@i OR a@i+1). Wrapped in
// `|-> 1'b0` (the same idiom first_match_seq.sv/seq_intersect.sv use)
// so the property is violated as soon as the composite sequence
// matches at all. `a`/`b`/`c` free: earliest match is b@1, c@2, a@1
// (or a@2), so the property should be violated at cycle 2 -- currently
// it is never reached: `within` is absent from the
// AssertionExprKind::Binary switches in
// assertion_expr_to_bool()/ltl_to_sat() (confirmed by inspection), so
// the default case returns a null Term and the whole assertion is
// silently dropped (propvec() stays empty) instead of being built or
// cleanly rejected.
module seq_within (input logic clk, input logic rst, input logic a, input logic b, input logic c);

  assert property (@(posedge clk) (a within (b ##1 c)) |-> 1'b0);

endmodule
