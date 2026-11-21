// `$cast` assigns its first argument the converted value of its second.
// A separate registration class in slang (CastTask), which is what
// decides whether the call reports having output arguments.
module system_task_cast (
    input logic       clk,
    input logic [7:0] a
);

  logic [7:0] x;

  initial x = 8'd0;

  always_ff @(posedge clk) begin
    $cast(x, a);
  end

  assert property (@(posedge clk) x == 8'd0);

endmodule
