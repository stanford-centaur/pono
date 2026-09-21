// A void task call used as a bare statement (`bump(a, b);`, no
// assignment) inside an always_ff block. The body is inlined and the
// output argument is written back to the caller's own variable,
// which is the whole point of a task.
//
// That copy-out is a blocking write inside a clocked block, so `b`
// is a flop: it holds the *previous* `a` plus one. The property says
// so, rather than the `b == a + 1` this fixture used to carry, which
// no run ever reached.
module void_task_call (
    input logic clk,
    input logic rst,
    input logic [3:0] a
);
  logic [3:0] b;

  task automatic bump(input logic [3:0] x, output logic [3:0] y);
    y = x + 4'd1;
  endtask

  always_ff @(posedge clk) begin
    if (rst) b <= 4'd0;
    else bump(a, b);
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> b == $past(a) + 4'd1);
endmodule
