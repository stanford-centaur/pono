// `initial forever @(posedge clk) ...` is a legacy structural way to
// spell what `always_ff @(posedge clk) ...` means: the timing-
// controlled body runs exactly once per clock edge, exactly like an
// always_ff block, so `q` should be a proper register (as opposed to a
// *bare* `forever` with no event control at all, which has no static
// iteration bound and is a genuine architectural boundary of this
// encoder's compile-time-unrolling model).
module forever_loop (input logic clk, input logic rst, input logic [3:0] din);

  logic [3:0] q;

  initial begin
    forever @(posedge clk) begin
      if (rst) q <= 4'd0;
      else q <= din;
    end
  end

  // q(T) is computed from rst/din as of cycle T-1, same as any
  // register.
  assert property (@(posedge clk) $past(rst) || q == $past(din));

endmodule
