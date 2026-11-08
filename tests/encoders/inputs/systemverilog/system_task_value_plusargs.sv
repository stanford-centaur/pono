// `$value$plusargs` writes the parsed plusarg into its second argument.
// A separate registration class in slang (MiscSystemFuncs), which is what
// decides whether the call reports having output arguments.
module system_task_value_plusargs (
    input logic       clk,
    input logic [7:0] a
);

  logic [7:0] x;

  initial x = 8'd0;

  always_ff @(posedge clk) begin
    $value$plusargs("v=%d", x);
  end

  assert property (@(posedge clk) x == 8'd0);

endmodule
