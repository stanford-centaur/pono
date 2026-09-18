// A plain user-defined SV `function` called with a symbolic (runtime-
// dependent) argument, distinct from a call with only compile-time-
// constant arguments (which slang's own constant evaluator resolves
// during elaboration, never reaching this encoder's expression walk).
// expr_to_term()'s Call case only recognizes a fixed list of system
// calls ($past/$stable/$changed/$rose/$fell/$onehot/$onehot0/
// $isunknown); any other call, including this one, throws
// "unsupported call", since user function bodies aren't inlined.
module user_function_call (input logic clk, input logic rst, input logic [3:0] a);

  function automatic logic [3:0] increment(logic [3:0] x);
    return x + 4'd1;
  endfunction

  logic [3:0] b;
  assign b = increment(a);

  assert property (@(posedge clk) b == a + 4'd1);

endmodule
