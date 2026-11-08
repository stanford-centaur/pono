// `$sformat` writes the formatted result into its first argument.
// A separate registration class in slang (StringFormatTask), which is what
// decides whether the call reports having output arguments.
module system_task_sformat (
    input logic       clk,
    input logic [7:0] a
);

  logic [7:0] x;
  string s;

  initial x = 8'd0;

  always_ff @(posedge clk) begin
    $sformat(s, "%0d", a);
  end

  assert property (@(posedge clk) x == 8'd0);

endmodule
