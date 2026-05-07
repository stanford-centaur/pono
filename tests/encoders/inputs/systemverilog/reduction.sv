// Tests reduction operators (&, |, ^).
// Property: !(all_ones && parity).  For a 3-bit input, all_ones=1
// requires a=3'b111, which has odd parity (parity=1).  BMC must
// falsify within 2 steps when a == 7.
module reduction (input logic clk, input logic [2:0] a);
  logic all_ones, parity, any_one;
  always_ff @(posedge clk) begin
    all_ones <= &a;
    any_one  <= |a;
    parity   <= ^a;
  end
  assert property (@(posedge clk) !(all_ones && parity));
endmodule
