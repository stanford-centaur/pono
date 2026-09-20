// weak() of `and`/`or` over multi-cycle sequences:
// offsets_ending_now() does not span these, so weak_seq_bool()
// cannot model them and the property falls through to the LTL
// tableau. Unwrapping it there hands each operand to
// try_strong_sequence(), which attaches the strong "must eventually
// complete" obligation -- the opposite of weak. It used to do that
// silently, making the assertion stronger than written; now the
// weak qualifier rides down the recursion and the sequence refuses.
module weak_seq_and_unsupported (
    input logic clk,
    input logic x,
    input logic y
);

  assert property (@(posedge clk) weak((x ##1 y) and (x ##1 y)));

endmodule
