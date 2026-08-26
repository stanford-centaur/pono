// Out-of-scope: a `specify` block (timing checks / path delays).
module specify_block (
    input  logic clk,
    input  logic rst,
    input  logic a,
    output logic y
);
  assign y = a;
  specify
    (a => y) = 2;
  endspecify
  assert property (@(posedge clk) 1'b1);
endmodule
