// An output argument makes the call write back to its caller, which
// an expression has no way to do.
module user_function_output_arg (
    input logic clk,
    input logic [3:0] a
);

  logic [3:0] o;
  logic r;

  function automatic logic f(input logic [3:0] x, output logic [3:0] y);
    y = x;
    return 1'b1;
  endfunction

  always_ff @(posedge clk) r <= f(a, o);

  assert property (@(posedge clk) r == r);

endmodule
