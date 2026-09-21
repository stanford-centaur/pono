// A concurrent assertion written inside a `program` block. The
// program's stimulus is simulation-only and stays unencoded -- a
// program drives the DUT along one particular scenario, and pinning
// the inputs to it would leave every other input sequence unexplored
// while still reporting a proof. An assertion is not stimulus,
// though, so it is encoded like any other property.
//
// This is the only assertion in the design, so if it were dropped the
// encoder would produce no properties at all rather than a passing
// one. It is false: the counter does reach 5.
program check (
    input logic clk,
    input logic [3:0] cnt
);

  assert property (@(posedge clk) cnt != 4'd5);

endprogram

module program_assertion (
    input logic clk,
    input logic rst
);

  logic [3:0] cnt;
  always_ff @(posedge clk) cnt <= rst ? 4'd0 : cnt + 4'd1;

  check c (.clk(clk), .cnt(cnt));

endmodule
