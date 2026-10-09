// Compound assignment (`|=`).  Slang expands `x |= y` into
// `x = LValueReference | y`, where LValueReference resolves to the
// current value of x.  Combining `|=` with a procedural for loop
// over the 4 bits of `din` builds an OR-reduction of din.
module compound_assign_or_reduce (input logic clk, input logic rst,
                                  input logic [3:0] din);
  logic any_set;
  logic reg_any_set;

  always_comb begin
    any_set = 1'b0;
    for (int i = 0; i < 4; ++i) begin
      any_set |= din[i];
    end
  end

  always_ff @(posedge clk) begin
    if (rst) reg_any_set <= 1'b0;
    else     reg_any_set <= any_set;
  end

  // Falsifiable: any non-zero din lets reg_any_set go high one
  // cycle after reset releases.  reg_any_set != 1 fails at cycle 2.
  assert property (@(posedge clk) reg_any_set != 1);
endmodule
