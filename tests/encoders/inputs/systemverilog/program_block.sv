// Out-of-scope: a verification `program` block.
program test_program (input logic clk);
  initial begin
    $display("hi");
  end
endprogram

module program_block (input logic clk, input logic rst);
  test_program tp (clk);
  assert property (@(posedge clk) 1'b1);
endmodule
