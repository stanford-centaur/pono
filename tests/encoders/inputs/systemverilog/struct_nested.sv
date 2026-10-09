// A struct nested inside another packed struct, with two independent
// fields of the inner struct updated at different rates so the test
// can pin an exact cycle where both cross specific thresholds
// simultaneously.  Exercises nested MemberAccess (`s.a.x`, `s.a.y`)
// on both sides of an always_ff assignment.
//
// x(k) = k-1, y(k) = 2*(k-1) for k >= 1, so x == 3 and y == 6 hold at
// exactly the same cycle, k = 4.
typedef struct packed {
  logic [3:0] x;
  logic [3:0] y;
} inner_t;

typedef struct packed {
  inner_t     a;
  logic [1:0] tag;
} outer_t;

module struct_nested (input logic clk, input logic rst);

  outer_t s;

  always_ff @(posedge clk) begin
    if (rst) begin
      s.a.x <= 4'd0;
      s.a.y <= 4'd0;
      s.tag <= 2'b00;
    end else begin
      s.a.x <= s.a.x + 4'd1;
      s.a.y <= s.a.y + 4'd2;
      s.tag <= s.tag + 2'd1;
    end
  end

  assert property (@(posedge clk) s.a.x != 4'd3 || s.a.y != 4'd6);

endmodule
