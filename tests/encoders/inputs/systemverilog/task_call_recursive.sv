// A task that calls itself has no point at which expansion stops.
module task_call_recursive (
    input logic clk,
    input logic [3:0] a
);

  logic [3:0] b;

  task automatic rec(input logic [3:0] x, output logic [3:0] y);
    rec(x, y);
  endtask

  always_ff @(posedge clk) rec(a, b);

  assert property (@(posedge clk) b == b);

endmodule
