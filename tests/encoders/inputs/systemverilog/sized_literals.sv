// Tests sized integer literals: hex, decimal, binary, and a negative
// signed literal that must be re-encoded as the unsigned bit pattern.
// Property: !(x == 8'hAB && y == 4'b1011 && z == -8'sd1).
// Each conjunct is independent, so BMC must falsify within 2 steps.
module sized_literals (input logic clk, input logic sel);
  logic [7:0] x;
  logic [3:0] y;
  logic signed [7:0] z;
  always_ff @(posedge clk) begin
    x <= sel ? 8'hAB : 8'h00;
    y <= sel ? 4'b1011 : 4'b0000;
    z <= sel ? -8'sd1 : 8'sd0;
  end
  assert property (@(posedge clk) !(x == 8'hAB && y == 4'b1011
                                    && z == -8'sd1));
endmodule
