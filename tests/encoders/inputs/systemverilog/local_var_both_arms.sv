// A subroutine local written under a runtime condition in *both*
// arms of an if/else. Between them the arms cover every path, so
// there is nothing undefined to merge against -- but the encoder
// binds one arm at a time and, finding no earlier binding to use for
// the other, refuses with "the local variable 'tmp' is written only
// under a runtime condition, so it has no value on the other path".
//
// It has a value on every path. What is missing is the merge of two
// arms that each supply one: the Ite at statement_encoder.cpp only
// knows how to combine a new write with a binding that already
// exists, so the first arm walked has nothing to pair with.
//
// The same refusal arrives through a function's own return variable
// (`if (s) pick = a; else pick = ...;`), which is the shape
// unbounded_recursion.sv trips on before its recursion is ever
// noticed -- so this bug is also what makes that one report the
// wrong thing.
//
// Intended: provable. `q = sel ? d : ~d` written without the call is.
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
