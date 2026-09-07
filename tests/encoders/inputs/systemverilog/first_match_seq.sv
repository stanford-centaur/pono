// `first_match`: restricts a (possibly multi-match) sequence to its
// earliest match.  `first_match(a ##[1:2] b) |-> 1'b0` is falsified as
// soon as the sequence matches at all (the consequent is always
// false), so this is really testing that `first_match` parses and
// composes with a range-delay sequence rather than changing whether
// at least one match exists.  Exact bound verified against the
// actual run rather than hand-derived.
module first_match_seq (input logic clk, input logic rst, input logic a, input logic b);

  assert property (@(posedge clk) first_match(a ##[1:2] b) |-> 1'b0);

endmodule
