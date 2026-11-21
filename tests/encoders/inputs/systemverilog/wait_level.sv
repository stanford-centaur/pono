// `wait (expr) statement`: level-sensitive event control. The
// condition is evaluated, and if it is false the statement stays
// blocked until the condition becomes true (LRM 9.4.3). Suspending
// part-way through a process needs a per-process resumption point,
// which a cycle-based transition system has no notion of.
//
// The grammar makes the guarded statement a *child* of the wait
// (Syntax 9-5), so this spelling and the `wait (go); x = ...;` one in
// wait_then.sv differ only in where the write sits in the tree, not
// in meaning: in both, `x` becomes 5 once `go` is high and not
// before.
module wait_level (
    input logic clk,
    input logic go
);

  logic [3:0] x;

  initial begin
    x = 4'd0;
    wait (go) x = 4'd5;
  end

  assert property (@(posedge clk) go || x == 4'd0);

endmodule
