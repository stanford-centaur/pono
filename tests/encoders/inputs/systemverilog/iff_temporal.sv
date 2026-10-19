// SVA `iff` with a genuinely temporal operand on one side (as
// opposed to two plain boolean operands, which route through a
// separate, already-safe fast path). `done` never holds, so
// `s_eventually done` is always false; `req` is tied high, so the
// two sides never match -- violated essentially immediately.
module iff_temporal (input logic clk, input logic rst);
  logic req, done;
  assign req = 1'b1;
  assign done = 1'b0;
  assert property (@(posedge clk) (s_eventually done) iff req);
endmodule
