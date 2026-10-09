// `expect (property_expr);`, a procedural blocking-wait statement
// distinct from `assert property (...)`. The ConcurrentAssertion
// statement handler explicitly recognizes AssertionKind::Expect and
// logs+skips it as a simulation-only construct.
module expect_property (input logic clk, input logic a);

  initial begin
    expect (@(posedge clk) a);
  end

endmodule
