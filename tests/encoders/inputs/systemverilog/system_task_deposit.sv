// `$deposit` sets a variable directly, outside any assignment.
// A separate registration class in slang (DepositTask), which is what
// decides whether the call reports having output arguments.
module system_task_deposit (
    input logic       clk,
    input logic [7:0] a
);

  logic [7:0] x;

  initial x = 8'd0;

  always_ff @(posedge clk) begin
    $deposit(x, 8'd9);
  end

  assert property (@(posedge clk) x == 8'd0);

endmodule
