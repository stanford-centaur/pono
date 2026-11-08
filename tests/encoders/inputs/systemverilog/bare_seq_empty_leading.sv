// weak_seq_empty_leading.sv without the qualifier `assert property`
// supplies for it anyway. Refuted at the same depth, or the two
// spellings of one property have come apart again.
module bare_seq_empty_leading (
    input logic clk,
    input logic x,
    input logic y
);

  assert property (@(posedge clk) x[*0:2] ##1 y);

endmodule
