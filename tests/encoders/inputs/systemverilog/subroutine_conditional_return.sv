// A `return` reached only under a runtime condition. The statements
// after it still run, so the value would be whichever assignment ran
// last -- which is why it is refused rather than guessed. Assigning
// the function's name under the same condition expresses the same
// thing and is modelled today, so this is a convenience rather than
// a missing semantics, but it is ordinary RTL style.
//
// Intended: provable, and the ternary spelling of the same function
// is.
module subroutine_conditional_return (input logic clk, input logic [7:0] d,
                                      output logic [7:0] q);

  function automatic [7:0] clamp(input [7:0] n);
    if (n > 8'd8) return 8'hFF;
    clamp = n;
  endfunction

  always_comb q = clamp(d);

  assert property (@(posedge clk) (d > 8'd8) ? q == 8'hFF : q == d);

endmodule
