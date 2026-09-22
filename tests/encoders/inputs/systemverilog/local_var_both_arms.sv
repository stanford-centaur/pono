// A subroutine local written under a runtime condition in *both*
// arms of an if/else, so between them the arms give it a value on
// every path and there is nothing undefined to merge against.
//
// The same shape arrives through a function's own return variable,
// `if (s) pick = a; else pick = ...;`.
//
// `q = sel ? d : ~d` written without the call holds, and so must
// this.
module local_var_both_arms (input logic clk, input logic [7:0] d,
                            input logic sel, output logic [7:0] q);

  function automatic [7:0] pick(input [7:0] a, input logic s);
    logic [7:0] tmp;
    if (s) tmp = a;
    else tmp = a ^ 8'hFF;
    pick = tmp;
  endfunction

  always_comb q = pick(d, sel);

  assert property (@(posedge clk) sel ? q == d : q == (d ^ 8'hFF));

endmodule
