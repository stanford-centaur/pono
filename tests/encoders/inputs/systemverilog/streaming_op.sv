// The streaming (bit-reversal) concatenation operator `{<<{...}}`.
// Reversing 4'b0001 yields 4'b1000; `a` is free, so BMC can pick
// a == 4'b0001 at the first post-reset transition, falsifying at
// cycle 2.
module streaming_op (input logic clk, input logic rst, input logic [3:0] a);

  logic [3:0] rev;

  always_ff @(posedge clk) begin
    if (rst) rev <= 4'd0;
    else rev <= {<<{a}};
  end

  assert property (@(posedge clk) rev != 4'b1000);

endmodule
