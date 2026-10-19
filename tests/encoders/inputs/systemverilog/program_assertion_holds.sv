// program_assertion.sv's twin with a property that does hold, so the
// refutation there is not just the encoder mangling a program-scope
// assertion into something unsatisfiable.
program check (
    input logic clk,
    input logic [3:0] cnt
);

  assert property (@(posedge clk) cnt <= 4'd15);

endprogram

module program_assertion_holds (
    input logic clk,
    input logic rst
);

  logic [3:0] cnt;
  always_ff @(posedge clk) cnt <= rst ? 4'd0 : cnt + 4'd1;

  check c (.clk(clk), .cnt(cnt));

endmodule
