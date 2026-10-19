// `a[*2:$]`. Ending now, a run of at least two ends with a run of
// exactly two, and where the run began changes nothing about
// whether it ends here -- so this is refuted at the same depth as
// `a[*2]` in unbounded_repeat_range_bounded.sv.
module unbounded_repeat_range (
    input logic clk,
    input logic rst,
    input logic a
);

  assert property (@(posedge clk) a[*2:$] |-> 1'b0);

endmodule
