// Several writes to one non-wire target in a single always_comb
// block. They are sequential: the later write lands on top of the
// earlier one. Constraining each write's own slice separately instead
// would assert p == 0 *and* p[3] == 1 at once, which is unsatisfiable
// and makes every property over the design pass vacuously --
// comb_partial_write_fails.sv is the paired refutation that rules
// that out.
//
// q's write is at a runtime index, which has no fixed slice to
// constrain at all and so has to compose the same way.
module comb_partial_write (
    input logic clk,
    input logic [2:0] i
);

  logic [7:0] p;
  logic [7:0] q;

  always_comb begin
    p = 8'd0;
    p[3] = 1'b1;
  end

  always_comb begin
    q = 8'd0;
    q[i] = 1'b1;
  end

  assert property (@(posedge clk) (p == 8'h08 && q == (8'd1 << i)));

endmodule
