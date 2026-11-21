// The same level-sensitive wait written with a null guarded statement
// and the work in the statements after it. "The procedural statements
// following the wait statement shall remain blocked until that
// condition becomes true" (LRM 9.4.3), so the write is as much
// governed by `go` here as in wait_level.sv.
module wait_then (
    input logic clk,
    input logic go
);

  logic [3:0] x;

  initial begin
    x = 4'd0;
    wait (go);
    x = 4'd5;
  end

  assert property (@(posedge clk) x == 4'd5);

endmodule
