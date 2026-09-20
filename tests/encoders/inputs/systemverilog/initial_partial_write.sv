// Two writes to one target in a single `initial` block, the second
// covering part of what the first did. They compose in order, so x
// starts at 8'h08.
//
// Constraining each write's own slice separately instead would
// assert x == 0 and x[3] == 1 at once, which no initial state
// satisfies, and an unsatisfiable initial state refutes nothing --
// initial_partial_write_fails.sv is the paired refutation that rules
// that out.
module initial_partial_write (
    input logic clk
);

  logic [7:0] x;

  initial begin
    x = 8'd0;
    x[3] = 1'b1;
  end

  always_ff @(posedge clk) x <= x;

  assert property (@(posedge clk) x == 8'h08);

endmodule
