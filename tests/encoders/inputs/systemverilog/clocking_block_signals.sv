// A clocking block's own variables, read (`cb.d`) and driven
// (`cb.o`). Each names the signal the block samples or drives rather
// than storage of its own, and with the default skews that happens at
// the clock edge -- the only moment this encoder has. So going
// through the block is going to the signal.
//
// The drive matters as much as the read: before `cb.o <= d` resolved
// to `o`, the write landed on the clocking variable instead and `o`
// was left undriven and free, with nothing said about it.
module clocking_block_signals (
    input logic clk,
    input logic [3:0] d
);

  logic [3:0] q, o;
  logic started;

  default clocking cb @(posedge clk);
    input d;
    output o;
  endclocking

  initial started = 1'b0;

  always_ff @(posedge clk) begin
    q <= cb.d;
    cb.o <= d;
    started <= 1'b1;
  end

  assert property (started |-> (q == $past(d)) and (o == $past(d)));

endmodule
