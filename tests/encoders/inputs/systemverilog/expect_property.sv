// `expect (property_expr);`, a procedural blocking-wait statement
// distinct from `assert property (...)` -- AssertionKind::Expect,
// which the ImmediateAssertion statement handler never checked for,
// silently dropping this construct entirely (no log, no throw, no
// encoding effect).
module expect_property (input logic clk, input logic a);

  initial begin
    expect (@(posedge clk) a);
  end

endmodule
