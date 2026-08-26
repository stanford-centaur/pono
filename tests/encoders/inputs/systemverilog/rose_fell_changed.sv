// `$rose`/`$fell`/`$changed`, checked as a genuine identity rather
// than one hand-picked trace: for a 1-bit signal, "changed" is
// exactly "rose or fell" at every cycle, for every possible history
// of `a` (a free input).
module rose_fell_changed (input logic clk, input logic rst, input logic a);

  assert property (@(posedge clk) ($rose(a) || $fell(a)) == $changed(a));

endmodule
