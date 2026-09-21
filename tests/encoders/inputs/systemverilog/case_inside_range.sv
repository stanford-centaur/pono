// `case (x) inside` matches by set membership, so its items may be
// ranges as well as values, and an integral item is compared with
// wildcard equality just as in the expression form -- `4'b1?01`
// leaves bit 2 free.
//
// Three selectors, each checked against the hand-written membership
// test one cycle later.
module case_inside_range (
    input logic clk,
    input logic rst,
    input logic [3:0] s
);

  logic yr, ym, yw;

  always_ff @(posedge clk) begin
    case (s) inside
      [4'd1:4'd3]: yr <= 1'b1;
      default: yr <= 1'b0;
    endcase
    case (s) inside
      4'd0, [4'd2:4'd3]: ym <= 1'b1;
      default: ym <= 1'b0;
    endcase
    case (s) inside
      4'b1?01: yw <= 1'b1;
      default: yw <= 1'b0;
    endcase
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (yr == $past(s >= 4'd1 && s <= 4'd3)
                        && ym == $past(s == 4'd0
                                       || (s >= 4'd2 && s <= 4'd3))
                        && yw == $past((s & 4'b1011) == 4'b1001)));

endmodule
