// Tests reduction operators (&, |, ^).
// Reset clears the registers; one cycle later, free input a = 3'b111
// produces all_ones=1 and parity=1 (3 ones is odd), falsifying the
// property.
module reduction (input logic clk, input logic rst, input logic [2:0] a);
  logic all_ones, parity, any_one;
  always_ff @(posedge clk) begin
    if (rst) begin
      all_ones <= 0;
      any_one  <= 0;
      parity   <= 0;
    end else begin
      all_ones <= &a;
      any_one  <= |a;
      parity   <= ^a;
    end
  end
  assert property (@(posedge clk) !(all_ones && parity));
endmodule
