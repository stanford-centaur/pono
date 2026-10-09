// Several calls in one expression, a call inside another call's
// argument, and a function that calls a second one. Each call binds
// the same formals afresh, so a binding left behind by one call
// would give another the wrong value.
module user_function_nested (
    input logic clk,
    input logic [3:0] a,
    input logic [3:0] b
);

  function automatic logic [3:0] inc(logic [3:0] x);
    return x + 4'd1;
  endfunction

  function automatic logic [3:0] inc2(logic [3:0] x);
    return inc(inc(x));
  endfunction

  assert property (@(posedge clk)
                   ((inc(a) + inc(b)) == (a + b + 4'd2)
                    && inc(inc(a)) == a + 4'd2
                    && inc2(b) == b + 4'd2));

endmodule
