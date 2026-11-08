// `$swrite` is `$write` aimed at a variable rather than the console.
// A separate registration class in slang (StringOutputTask), which is what
// decides whether the call reports having output arguments.
module system_task_swrite (
    input logic       clk,
    input logic [7:0] a
);

  logic [7:0] x;
  string s;

  initial x = 8'd0;

  always_ff @(posedge clk) begin
    $swrite(s, "%0d", a);
  end

  assert property (@(posedge clk) x == 8'd0);

endmodule
