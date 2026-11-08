// `$fscanf` parses into its output arguments. It is a function rather
// than a task, so a design is free to discard its result, which leaves
// those writes as the only thing the call was for.
// A separate registration class in slang (NonConstFuncs), which is what
// decides whether the call reports having output arguments.
module system_task_fscanf (
    input logic       clk,
    input logic [7:0] a
);

  logic [7:0] x;
  integer fd;

  initial x = 8'd0;

  always_ff @(posedge clk) begin
    $fscanf(fd, "%d", x);
  end

  assert property (@(posedge clk) x == 8'd0);

endmodule
