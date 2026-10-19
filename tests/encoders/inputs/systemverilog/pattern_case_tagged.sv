// A `tagged` pattern matches on a tagged union's discriminant, which
// a *packed* one does have at a defined position: LRM 7.3.2 puts the
// tag in the top bits, holding the member's declaration index, and
// right-justifies each member below it. The bits between are
// undefined, and nothing reads them.
//
// Two members here, so one tag bit: `raw[8]` picks the arm, and both
// members read from the low bits regardless of which is live. The
// twin below refutes the reading where the tag is ignored and the
// payload alone decides.
typedef union tagged packed {
  logic [3:0] A;
  logic [7:0] B;
} pattern_case_tagged_u;

module pattern_case_tagged (
    input logic clk,
    input logic [8:0] raw
);

  pattern_case_tagged_u u;
  logic [7:0] y;

  always_comb begin
    u = raw;
    case (u) matches
      tagged A .n: y = {4'd0, n};
      tagged B .n: y = n;
      default: y = 8'd0;
    endcase
  end

  assert property (@(posedge clk)
                   ((raw[8] == 1'b0) |-> (y == {4'd0, raw[3:0]}))
                   and ((raw[8] == 1'b1) |-> (y == raw[7:0])));

endmodule
