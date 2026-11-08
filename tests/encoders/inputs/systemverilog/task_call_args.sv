// The argument directions a task write-back has to get right: two
// separate outputs, an inout that is read and then written, and a
// call made only on one path, whose write-back must be guarded by
// that same condition rather than happening unconditionally.
module task_call_args (
    input logic clk,
    input logic rst,
    input logic en,
    input logic [3:0] a
);

  logic [3:0] lo, hi, c, g;

  task automatic split(input logic [3:0] x, output logic [3:0] p,
                       output logic [3:0] q);
    p = x + 4'd1;
    q = x + 4'd2;
  endtask

  task automatic step(inout logic [3:0] v);
    v = v + 4'd1;
  endtask

  task automatic bump(input logic [3:0] x, output logic [3:0] y);
    y = x + 4'd1;
  endtask

  always_ff @(posedge clk) begin
    if (rst) begin
      lo <= 4'd0;
      hi <= 4'd0;
      c  <= 4'd0;
      g  <= 4'd0;
    end else begin
      split(a, lo, hi);
      step(c);
      if (en) bump(a, g);
    end
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (lo == $past(a) + 4'd1
                        && hi == $past(a) + 4'd2
                        && c == $past(c) + 4'd1
                        && g == ($past(en) ? $past(a) + 4'd1 : $past(g))));

endmodule
