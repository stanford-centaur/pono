// A void task call used as a bare statement (`bump(a, b);`, no
// assignment) inside an always_ff block. Task bodies aren't inlined
// by this encoder (the same limitation Gap_UserFunctionCall documents
// for a plain function call in expression context), so the task's
// side effect on `b` is not applied -- process_statement() throws
// "unsupported call" rather than silently dropping it, so the
// fixture's own `b == a + 1` invariant can't be checked either way.
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

  assert property (@(posedge clk) rst || b == a + 4'd1);
endmodule
