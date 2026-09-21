// A function that calls itself: each expansion binds the same
// formals again, so there is no point at which it stops.
module user_function_recursive (
    input logic clk,
    input logic [3:0] a
);

  function automatic logic [3:0] rec(logic [3:0] x);
    rec = rec(x);
  endfunction

  assert property (@(posedge clk) rec(a) == a);

endmodule
