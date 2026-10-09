// `while` as a compile-time-bounded loop, mirroring the OR-reduction
// that the already-supported `for` loop performs in compound_assign.sv
// -- same intent, different loop keyword.  `din` is free, so BMC can
// set any bit and falsify at cycle 2, *if* the while loop actually
// runs.
module while_loop (input logic clk, input logic rst, input logic [3:0] din);

  logic any_set, reg_any_set;

  always_comb begin
    int i;
    any_set = 1'b0;
    i = 0;
    while (i < 4) begin
      any_set |= din[i];
      i = i + 1;
    end
  end

  always_ff @(posedge clk) begin
    if (rst) reg_any_set <= 1'b0;
    else reg_any_set <= any_set;
  end

  assert property (@(posedge clk) reg_any_set != 1'b1);

endmodule
