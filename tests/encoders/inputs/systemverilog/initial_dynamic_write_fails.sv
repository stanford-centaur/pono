// The other half of initial_dynamic_write.sv: only the selected bit
// is pinned, so claiming the whole word is that one-hot value must
// be refuted. Splicing onto a zero seed rather than the variable's
// own term would have made this hold.
module initial_dynamic_write_fails (
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

  assert property (@(posedge clk) (!seen) |-> r == (8'd1 << i));

endmodule
