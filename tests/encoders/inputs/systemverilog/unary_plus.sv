// Unary `+` (a no-op per the LRM: `+a` == `a`), registered.  `a` is
// free, so BMC can pick a == 5 at the first post-reset transition,
// falsifying at cycle 2.
module unary_plus (input logic clk, input logic rst, input logic [3:0] a);

  logic [3:0] p;

  always_ff @(posedge clk) begin
    if (rst) p <= 4'd0;
    else p <= +a;
  end

  assert property (@(posedge clk) p != 4'd5);

endmodule
