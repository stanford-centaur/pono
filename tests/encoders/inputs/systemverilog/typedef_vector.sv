// A `typedef` over a plain bit-vector (rather than a struct) used for
// both a state variable's type and its increment step size, wired
// through a real register over several cycles.
//
// val(k) = 5*(k-1) for k >= 1, so val == 20 first at cycle 5.
typedef logic [7:0] byte_t;
typedef byte_t step_t;

module typedef_vector (input logic clk, input logic rst);

  byte_t val;
  step_t step;

  assign step = 8'd5;

  always_ff @(posedge clk) begin
    if (rst) val <= 8'd0;
    else val <= val + step;
  end

  assert property (@(posedge clk) val != 8'd20);

endmodule
