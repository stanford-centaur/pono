// The packed qualifier is what makes a tagged union modelable: LRM
// 7.3 gives an unpacked union "no required representation" at all, so
// its tag has no position to be read from and no layout to match
// against. Refused at the type, before any pattern is reached.
typedef union tagged {
  void Invalid;
  int  Valid;
} unpacked_tagged_t;

module unpacked_tagged_union (
    input logic clk
);

  unpacked_tagged_t u;
  logic [3:0] y;

  always_comb begin
    case (u) matches
      tagged Valid .n: y = 4'd1;
      tagged Invalid: y = 4'd0;
      default: y = 4'd2;
    endcase
  end

  assert property (@(posedge clk) 1'b1);

endmodule
