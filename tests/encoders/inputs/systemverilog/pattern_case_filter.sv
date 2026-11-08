// A `&&&` filter narrows an item beyond its pattern, and may read
// the names that pattern just bound -- so the binding has to be in
// scope while the filter is built, not only while the arm's
// statement is walked. With the filter present the item is no
// longer irrefutable, which is why this one needs its `default`.
module pattern_case_filter (input logic clk, input logic [3:0] x);
  logic [3:0] y;
  always_comb begin
    case (x) matches
      .v &&& (v > 4'd7) : y = 4'd1;
      default           : y = 4'd0;
    endcase
  end
  assert property (@(posedge clk) (x > 4'd7) == (y == 4'd1));
endmodule
