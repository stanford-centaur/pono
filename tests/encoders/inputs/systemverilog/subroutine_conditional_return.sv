// A `return` reached only under a runtime condition, after which the
// remaining statements still run -- so the value is whichever
// assignment ran last. Assigning the function's name under the same
// condition says the same thing, but this is ordinary RTL style.
//
// The ternary spelling of the same function holds, and so must this.
module subroutine_conditional_return (input logic clk, input logic [7:0] d,
                                      output logic [7:0] q);

  function automatic [7:0] clamp(input [7:0] n);
    if (n > 8'd8) return 8'hFF;
    clamp = n;
  endfunction

  always_comb q = clamp(d);

  assert property (@(posedge clk) (d > 8'd8) ? q == 8'hFF : q == d);

endmodule
