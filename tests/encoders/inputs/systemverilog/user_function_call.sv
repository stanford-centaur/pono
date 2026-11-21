// A plain user-defined SV `function` called with a symbolic (runtime-
// dependent) argument, distinct from a call with only compile-time-
// constant arguments, which slang's own constant evaluator resolves
// during elaboration and which never reaches this encoder's
// expression walk at all. The body is inlined where the call
// appears, so the assertion below pins the value it produces rather
// than just that it produces one.
module user_function_call (input logic clk, input logic rst, input logic [3:0] a);

  function automatic logic [3:0] increment(logic [3:0] x);
    return x + 4'd1;
  endfunction

  logic [3:0] b;
  assign b = increment(a);

  assert property (@(posedge clk) b == a + 4'd1);

endmodule
