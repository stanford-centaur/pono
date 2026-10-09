// s_always_range.sv with the window widened by one cycle: four cycles
// after a 0 the wrapping counter is back at 0, so the property fails.
// The two fixtures differ only in the window's upper bound, so this
// pair is what pins down that the bound is actually honoured.
module s_always_range_fails (
    input logic clk,
    input logic rst
);

  logic [1:0] count;

  always_ff @(posedge clk) begin
    if (rst) count <= 2'd0;
    else count <= count + 2'd1;
  end

  assert property (@(posedge clk)
                    (!rst && count == 2'd0) |-> s_always [1:4] (count != 2'd0));

endmodule
