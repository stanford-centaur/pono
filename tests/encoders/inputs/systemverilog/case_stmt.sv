// Tests a case statement inside always_ff.
// When s == 3, x is set to 8.
// Property: x != 8.  BMC must falsify within 2 steps.
module case_stmt (input logic clk, input logic [1:0] s);
  logic [3:0] x;
  always_ff @(posedge clk) begin
    case (s)
      0: x <= 1;
      1: x <= 2;
      2: x <= 4;
      3: x <= 8;
    endcase
  end
  assert property (@(posedge clk) x != 8);
endmodule
