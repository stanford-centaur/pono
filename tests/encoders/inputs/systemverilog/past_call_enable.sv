// `$past(expr, n, enable)`'s `enable` argument gates whether a given
// cycle's sample is taken at all: it should freeze the whole delayed
// history on a disabled cycle, not just skip that cycle's value.
// Checked as a genuine identity for every possible `data`/`en` history
// (like rose_fell_changed.sv), not one hand-picked trace: `held` mirrors
// $past's own enable-gated freeze/update rule directly (sample `data`
// when `en` holds, otherwise hold), and both start from the same
// deterministic zero state (`held` via `initial`, $past's internal
// latch is always zero-initialized), so a correct $past(data, 1, en)
// must agree with `held` at every cycle, with no reset/warm-up window
// needed.
module past_call_enable (input logic clk, input logic en,
                         input logic [3:0] data);
  logic [3:0] held;
  initial held = '0;
  always_ff @(posedge clk)
    if (en) held <= data;

  assert property (@(posedge clk) $past(data, 1, en) == held);
endmodule
