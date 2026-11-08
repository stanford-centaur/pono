// A `case` whose arms cover every value of the selector, with no
// `default`. The definite-assignment scan is syntactic and reads a
// missing `default` as leaving a path open, so it marks `y` a latch;
// believing that would make `y` a register and delay it a cycle.
//
// Asking the solver whether the accumulated value can still be
// influenced by `y`'s own previous value settles it: `s` is one bit
// and both values are matched, so the fallback is unreachable and `y`
// is plain combinational. The delayed twin of this fixture is what
// stops the register reading coming back.
module case_full_coverage (
    input logic clk,
    input logic s,
    input logic [3:0] p,
    input logic [3:0] q
);

  logic [3:0] y;

  always_comb begin
    case (s)
      1'b0: y = p;
      1'b1: y = q;
    endcase
  end

  assert property (@(posedge clk) y == (s ? q : p));

endmodule
