// A level-sensitive latch: `q` should hold its previous value whenever
// `en` (and `rst`) are both low, and load `d` whenever `en` is high.
// Modeled the same way a register is (a state variable whose next
// value defaults to holding itself), just triggered by blocking `=`
// inside always_latch instead of nonblocking `<=` inside always_ff --
// this encoder's per-cycle abstraction has no notion of "edge" vs.
// "level" triggering, so the two should be indistinguishable here.
module always_latch_test (
    input logic clk,
    input logic rst,
    input logic en,
    input logic [3:0] d
);

  logic [3:0] q;

  always_latch begin
    if (rst) q = 4'd0;
    else if (en) q = d;
  end

  assert property (@(posedge clk) (rst || en) || (q == $past(q)));

endmodule
