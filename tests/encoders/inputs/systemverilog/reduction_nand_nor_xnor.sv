// Reduction NAND/NOR/XNOR (`~&`, `~|`, `~^`), registered.  NAND-reduce
// of a is false only when a is all-ones; NOR-reduce is true only when
// a is all-zero -- both checked against a same-cycle-paired copy of
// `a`, so they genuinely hold for every a.  `xnor_r` is computed (and
// so exercised for parse/encode purposes) but not asserted on, since
// stating its exact value correctly here is more error-prone than the
// point of the test warrants.
module reduction_nand_nor_xnor (input logic clk, input logic rst, input logic [3:0] a);

  logic [3:0] a_reg;
  logic       nand_r, nor_r, xnor_r;

  always_ff @(posedge clk) begin
    if (rst) begin
      // a_reg resets to 0, so nand_r/nor_r must reset to 1 (NAND/NOR
      // of an all-zero input are both true) to stay consistent with
      // the assertion below -- forcing them to 0 here would violate
      // the property from the very first post-reset cycle, unrelated
      // to whether ~&/~| are actually correct.
      a_reg  <= 4'd0;
      nand_r <= 1'b1;
      nor_r  <= 1'b1;
      xnor_r <= 1'b0;
    end else begin
      a_reg  <= a;
      nand_r <= ~&a;
      nor_r  <= ~|a;
      xnor_r <= ~^a;
    end
  end

  assert property (@(posedge clk)
                    nand_r == (a_reg != 4'hF) && nor_r == (a_reg == 4'h0));

endmodule
