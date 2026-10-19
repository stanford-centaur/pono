// A `forever` loop with no event control at all -- unlike
// `initial forever @(posedge clk) ...` (a legacy structural spelling
// of `always @(posedge clk) ...`, which this encoder does support),
// this has no static iteration bound whatsoever and can't be unrolled
// by the compile-time-bounded model. Must throw a clear error.
module bare_forever (input logic clk, input logic rst, input logic [3:0] din);

  logic [3:0] q;

  initial begin
    forever begin
      q = din;
    end
  end

  assert property (@(posedge clk) rst || q == din);

endmodule
