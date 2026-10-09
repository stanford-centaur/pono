// The twin of clocking_block_signals.sv, claiming the registers keep
// up with `d` in the same cycle. They do not, and this has to be
// refuted -- otherwise `cb.d` resolving to something free rather than
// to `d` would go unnoticed, since an unconstrained value satisfies
// nothing in particular but fails nothing in particular either.
module clocking_block_signals_delayed (
    input logic clk,
    input logic [3:0] d
);

  logic [3:0] q;
  logic started;

  default clocking cb @(posedge clk);
    input d;
  endclocking

  initial started = 1'b0;

  always_ff @(posedge clk) begin
    q <= cb.d;
    started <= 1'b1;
  end

  assert property (started |-> q == d);

endmodule
