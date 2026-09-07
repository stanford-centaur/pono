// `a[*2:$]` (an explicit unbounded consecutive-repetition range), same
// architectural boundary as unbounded_repeat_star.sv -- distinct
// syntax, same `[*n:$]` unbounded-max shape.
module unbounded_repeat_range (input logic clk, input logic rst, input logic a);

  assert property (@(posedge clk) a[*2:$] |-> 1'b0);

endmodule
