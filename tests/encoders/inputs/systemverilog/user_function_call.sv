// A plain user-defined SV `function` called with a symbolic (runtime-
// dependent) argument -- a common RTL idiom (a small combinational
// helper function), distinct from a function called only with
// compile-time-constant arguments (which slang's own constant
// evaluator resolves during elaboration, e.g. in a parameter
// expression, never reaching this encoder's own expression walk at
// all). expr_to_term()'s Call case only recognizes a fixed list of
// system calls ($past/$stable/$changed/$rose/$fell/$onehot/$onehot0/
// $isunknown); any other call -- including this one -- throws
// "unsupported call", since inlining a user function's body isn't
// implemented.
module user_function_call (input logic clk, input logic rst, input logic [3:0] a);

  function automatic logic [3:0] increment(logic [3:0] x);
    return x + 4'd1;
  endfunction

  logic [3:0] b;
  assign b = increment(a);

  assert property (@(posedge clk) b == a + 4'd1);

endmodule
