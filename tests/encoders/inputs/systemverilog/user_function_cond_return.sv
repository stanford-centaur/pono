// A `return` reached only under a runtime condition. The statements
// after it are still walked, so the value would be whatever the last
// reached assignment left rather than this one -- assigning the
// function's name under the condition says the same thing and is
// modelled, so this is refused rather than silently mis-encoded.
module user_function_cond_return (
    input logic clk,
    input logic [3:0] a
);

  function automatic logic [3:0] pick(logic [3:0] x);
    if (x[0]) return 4'd1;
    return 4'd2;
  endfunction

  assert property (@(posedge clk) pick(a) == 4'd1);

endmodule
