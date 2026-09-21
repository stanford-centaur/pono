// A `&&&`-joined multi-condition `if` (`if (a &&& b) ...`, LRM 12.4.4)
// must AND together every joined condition, not just the first --
// mirrors multi_cond_ternary.sv's ConditionalOp coverage for the
// ConditionalStatement path.
module multi_cond_if (input logic clk,
                      input logic a,
                      input logic b);

  logic [3:0] result;

  always_comb begin
    if (a &&& b) result = 4'd1;
    else result = 4'd0;
  end

  assert property (@(posedge clk) !(a && !b) || (result == 4'd0));

endmodule
