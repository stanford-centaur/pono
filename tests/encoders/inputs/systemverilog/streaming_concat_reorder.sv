// The `<<` half of streaming concatenation, in both directions, with
// a slice size that does *not* divide the width -- the case where
// packing and unpacking are different permutations rather than the
// same one applied twice.
//
// Per LRM 11.4.14.2 the blocks are cut from the right-hand end, so a
// 6-bit value under `<<4` is a 2-bit left-over `L` above one 4-bit
// block `B`; reversing puts `B` first and leaves `L` at the bottom.
// Unpacking (11.4.14.3) has to run that backwards, cutting the short
// piece off the *bottom* and returning it to the top.
//
// The third target streams two expressions rather than one, which a
// single-operand encoding could not express at all.
module streaming_concat_reorder (
    input logic clk,
    input logic [5:0] a,
    input logic [3:0] p,
    input logic [3:0] q
);

  logic [5:0] packed_ragged;
  logic [5:0] unpacked_ragged;
  logic [7:0] two_streams;
  logic started;

  initial started = 1'b0;

  always_ff @(posedge clk) begin
    // Pack: a == L ++ B, so the stream is B ++ L.
    packed_ragged <= {<<4{a}};
    // Unpack: the inverse, so the target is a[1:0] ++ a[5:2].
    {<<4{unpacked_ragged}} <= a;
    two_streams <= {<<4{p, q}};
    started <= 1'b1;
  end

  assert property (@(posedge clk) started
                   |-> (packed_ragged == {$past(a[3:0]), $past(a[5:4])}
                        && unpacked_ragged[5:4] == $past(a[1:0])
                        && unpacked_ragged[3:0] == $past(a[5:2])
                        && two_streams == {$past(q), $past(p)}));

endmodule
