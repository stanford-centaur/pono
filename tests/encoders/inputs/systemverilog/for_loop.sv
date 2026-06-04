// Procedural `for` loop in an always_comb block.  The loop
// computes the popcount of a 4-bit input by summing into a scalar
// (NamedValue LHS, since bit-select LHS is not yet supported).
//
// After reset, popcount[t] == popcount(din[t-1]) <= 4 at every
// cycle, so the property `popcount != 5` first fails at cycle... it
// can't actually fail, popcount is bounded by 4.  Make the test
// falsifiable by asserting a tighter bound: popcount <= 3.
module for_loop_popcount (input logic clk, input logic rst,
                          input logic [3:0] din);

  logic [3:0] popcount;
  logic [3:0] reg_popcount;

  always_comb begin
    popcount = '0;
    for (int i = 0; i < 4; ++i) begin
      popcount = popcount + din[i];
    end
  end

  always_ff @(posedge clk) begin
    if (rst) reg_popcount <= '0;
    else     reg_popcount <= popcount;
  end

  // Falsifiable: din == 4'b1111 forces popcount == 4, so reg_popcount
  // reaches 4 one cycle after reset releases.
  assert property (@(posedge clk) reg_popcount != 4);

endmodule
