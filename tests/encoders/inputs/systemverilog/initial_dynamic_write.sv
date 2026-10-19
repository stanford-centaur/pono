// A runtime-indexed write in an `initial` block. It needs no fixed
// slice after all: the splice composes onto the variable's own
// term, so the constraint reads "r equals r with these bits
// replaced" -- a tautology everywhere else and a definition exactly
// where the write lands. That pins the selected bit and leaves the
// rest of the initial value free.
//
// `i` is a state variable, so the design knows it at time 0; an
// initial value may not depend on an input, which has no value
// then. `seen` is initialised in the same block so that `!seen`
// really does single out time 0 rather than being free there.
//
// initial_dynamic_write_fails.sv pins the "rest stays free" half.
module initial_dynamic_write (
    input logic clk
);

  logic [7:0] r;
  logic [2:0] i;
  logic seen;

  always_ff @(posedge clk) begin
    i <= i;
    seen <= 1'b1;
  end

  initial begin
    seen = 1'b0;
    r[i] = 1'b1;
  end

  assert property (@(posedge clk) (!seen) |-> r[i] == 1'b1);

endmodule
