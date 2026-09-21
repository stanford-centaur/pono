// Building a tagged union value (`tagged B (v)`), the other half of
// matching one. The same LRM 7.3.2 layout in the other direction: the
// member's declaration index in the top bits, the value
// right-justified below, and the gap between them left unconstrained
// rather than zeroed, since the standard gives those bits no value.
//
// `hit_b` checks the round trip, and `hit_other` checks that a value
// built as `A` does not match the `B` arm -- without the tag, the
// payload bits alone would let it.
typedef union tagged packed {
  logic [3:0] A;
  logic [7:0] B;
} tagged_union_construct_u;

module tagged_union_construct (
    input logic clk,
    input logic [7:0] v,
    input logic [3:0] w
);

  tagged_union_construct_u ub, ua;
  logic [7:0] y;
  logic hit_b, hit_other;

  always_comb begin
    ub = tagged B (v);
    ua = tagged A (w);
    y = 8'd0;
    hit_b = 1'b0;
    hit_other = 1'b0;
    case (ub) matches
      tagged B .n: begin
        y = n;
        hit_b = 1'b1;
      end
      default: ;
    endcase
    case (ua) matches
      tagged B .n: hit_other = 1'b1;
      default: ;
    endcase
  end

  assert property (@(posedge clk) hit_b and (y == v) and (!hit_other));

endmodule
