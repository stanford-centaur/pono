// Out-of-scope: a DPI-C imported function called from RTL.
module dpi_import (input logic clk, input logic rst, input int a);
  import "DPI-C" function int c_add(int x, int y);
  int result;
  always_ff @(posedge clk) begin
    if (rst) result <= 0;
    else result <= c_add(a, 1);
  end
  assert property (@(posedge clk) 1'b1);
endmodule
