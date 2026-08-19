// `initial forever @(posedge clk) ...` is a legacy structural way to
// spell what `always_ff @(posedge clk) ...` means -- but `forever` has
// no static iteration bound at all (unlike `for`/`while`/`repeat`,
// which this encoder unrolls up to a compile-time-computable count), a
// genuine architectural boundary of the compile-time-unrolling model,
// not a "not implemented yet" gap. Must throw a clear error rather
// than silently dropping the assignment inside it.
module forever_loop (input logic clk, input logic rst, input logic [3:0] din);

  logic [3:0] q;

  initial begin
    forever @(posedge clk) begin
      if (rst) q <= 4'd0;
      else q <= din;
    end
  end

  assert property (@(posedge clk) rst || q == din);

endmodule
