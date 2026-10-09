// Out-of-scope: a dynamic array (`int arr[];`), unlike the
// fixed-size packed arrays covered in test_systemverilog_types.cpp.
module dynamic_array (input logic clk, input logic rst);
  int arr[];
  always_ff @(posedge clk) begin
    if (rst) arr = new[4];
  end
  assert property (@(posedge clk) 1'b1);
endmodule
