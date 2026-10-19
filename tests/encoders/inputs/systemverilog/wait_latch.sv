// The level-sensitive latch idiom: a free-running `always` whose body
// is a wait. `q` follows `d` while `en` is high and holds otherwise.
// Expressing it needs the process to park at the wait between passes,
// which is the per-process state a transition system does not carry.
module wait_latch (
    input logic clk,
    input logic en,
    input logic [3:0] d
);

  logic [3:0] q;

  initial q = 4'd0;

  always begin
    wait (en) q = d;
  end

  assert property (@(posedge clk) en || q == $past(q));

endmodule
