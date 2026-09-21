// A runtime-indexed write inside a child whose output port is
// connected to a concatenation, so the port's bits are spread across
// two parent-side signals. Which of them the write reaches is only
// known at runtime.
//
// Index 0..3 lands in `lo` and 4..7 in `hi`, and the half the write
// misses must keep what it had -- a guard that were too generous
// would disturb both. dynamic_write_port_alias_fails.sv is the
// paired refutation, since asserting only that the right half
// changed would also pass if neither did.
module dynamic_write_port_alias (
    input logic clk,
    input logic rst,
    input logic [2:0] i
);

  logic [3:0] hi, lo;

  alias_driver u (
      .clk(clk),
      .rst(rst),
      .i  (i),
      .p  ({hi, lo})
  );

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> ($past(i) < 3'd4
                        ? (lo[$past(i)] == 1'b1 && hi == $past(hi))
                        : (hi[$past(i)-3'd4] == 1'b1 && lo == $past(lo))));

endmodule

module alias_driver (
    input logic clk,
    input logic rst,
    input logic [2:0] i,
    output logic [7:0] p
);
  always_ff @(posedge clk) begin
    if (rst) p <= 8'd0;
    else p[i] <= 1'b1;
  end
endmodule
