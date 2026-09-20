// `a[*]` (`[*0:$]`, unbounded consecutive repetition) used as the
// antecedent of `|->`, the code path that handles bounded sequences
// generally (see offsets_ending_now() in frontends/systemverilog/
// sva.cpp). `a[*]` matches emptily, and on its own that is all the
// antecedent needs to match -- an empty match occupies no cycles, so
// there is no cycle at which the antecedent could be said to end and
// hand over to the consequent. Must throw rather than silently
// dropping the assertion, or reading it as an obligation on every
// cycle that the `a` would then play no part in.
module unbounded_repeat_star (input logic clk, input logic rst, input logic a);

  assert property (@(posedge clk) a[*] |-> 1'b0);

endmodule
