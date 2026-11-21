// The two idiomatic ways to clear a memory on reset, in one design:
// `foreach`, which unrolls to one constant-index element write per
// entry, and a whole-array constant pattern, which becomes a single
// constant array rather than a chain of writes.
//
// Also covers a non-zero-based declared range: element 0 of the SMT
// array is `mem[3]`, so every read and write has to subtract the
// declared lower bound.  An off-by-one there would break the
// read-after-write invariant below.
module unpacked_array_reset (
    input logic clk,
    input logic rst,
    input logic [7:0] wdata
);

  logic [7:0] mem[3:18];
  logic [7:0] shadow[0:3];

  always_ff @(posedge clk) begin
    if (rst) begin
      foreach (mem[i]) mem[i] <= 8'd0;
      shadow <= '{default: 0};
    end else begin
      mem[3] <= wdata;
      shadow[0] <= wdata;
    end
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> mem[3] == $past(wdata) && shadow[0] == $past(wdata));

endmodule
