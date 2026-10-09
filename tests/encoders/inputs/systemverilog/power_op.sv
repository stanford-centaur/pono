// The `**` power operator, registered.  `a` is a free 3-bit input
// (0..7), so a**2 ranges 0..49; BMC can pick a==3 at the first
// post-reset transition (a**2==9), falsifying at cycle 2.
module power_op (input logic clk, input logic rst, input logic [2:0] a);

  logic [7:0] result;

  always_ff @(posedge clk) begin
    if (rst) result <= 8'd0;
    else result <= a ** 2;
  end

  assert property (@(posedge clk) result != 8'd9);

endmodule
