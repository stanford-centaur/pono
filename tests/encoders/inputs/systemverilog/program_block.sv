// A `program` whose body is stimulus only: a `$display` in an
// `initial`, with no per-cycle meaning. Ignored rather than
// encoded or rejected. A program containing an assertion is a
// different matter -- see program_assertion.sv.
program test_program (input logic clk);
  initial begin
    $display("hi");
  end
endprogram

module program_block (input logic clk, input logic rst);
  test_program tp (clk);
  assert property (@(posedge clk) 1'b1);
endmodule
