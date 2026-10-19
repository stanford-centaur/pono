// `break` guarded by a runtime signal (`cond`, a free input) rather
// than a compile-time-constant loop-counter comparison. The
// Conditional case's constant-fold fast path only applies when the
// guard is const-evaluable, so this falls back to the general
// symbolic-guard path, which processes both arms of every `if`
// unconditionally -- meaning the `break` is always reached while
// walking the design, with no way to correctly model "only break when
// cond is true" as C++-level control flow. This must throw a clear
// error rather than silently doing nothing or (worse) always/never
// breaking regardless of `cond`.
module break_runtime_dependent (input logic clk, input logic rst, input logic cond);

  logic [3:0] din;
  logic       any_set, reg_any_set;

  assign din = 4'b1000;

  always_comb begin
    any_set = 1'b0;
    for (int i = 0; i < 4; ++i) begin
      if (cond) break;
      any_set |= din[i];
    end
  end

  always_ff @(posedge clk) begin
    if (rst) reg_any_set <= 1'b0;
    else reg_any_set <= any_set;
  end

  assert property (@(posedge clk) reg_any_set != 1'b1);

endmodule
