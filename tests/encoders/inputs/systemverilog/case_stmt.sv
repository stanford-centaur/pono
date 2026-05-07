// Tests a case statement inside always_ff.
// Reset clears x; one cycle later, s == 3 selects the case that sets
// x to 8, falsifying the property.
module case_stmt (input logic clk, input logic rst, input logic [1:0] s);
  logic [3:0] x;
  always_ff @(posedge clk) begin
    if (rst) x <= 0;
    else case (s)
      0: x <= 1;
      1: x <= 2;
      2: x <= 4;
      3: x <= 8;
    endcase
  end
  assert property (@(posedge clk) x != 8);
endmodule
