// `a[*]` (`[*0:$]`, unbounded consecutive repetition) used as the
// antecedent of `|->`, the code path that handles bounded sequences
// generally (see offsets_ending_now() in frontends/systemverilog/
// sva.cpp). An unbounded repeat count can't be unrolled by this
// encoder's compile-time-bounded model -- must throw a clear error
// rather than silently dropping the assertion.
module unbounded_repeat_star (input logic clk, input logic rst, input logic a);

  assert property (@(posedge clk) a[*] |-> 1'b0);

endmodule
