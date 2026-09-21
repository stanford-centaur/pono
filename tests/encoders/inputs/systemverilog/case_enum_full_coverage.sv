// The shape this actually shows up as: an enum selector with every
// state enumerated and no `default`, decoded in an always_comb. It is
// ordinary RTL rather than a corner, which is why the syntactic
// reading of "no default" mattered so much.
typedef enum logic [1:0] {
  IDLE,
  RUN,
  DONE,
  HALT
} case_enum_st_t;

module case_enum_full_coverage (
    input logic clk,
    input case_enum_st_t st
);

  logic [3:0] y;

  always_comb begin
    case (st)
      IDLE: y = 4'd0;
      RUN:  y = 4'd1;
      DONE: y = 4'd2;
      HALT: y = 4'd3;
    endcase
  end

  assert property (@(posedge clk) (st == IDLE) |-> (y == 4'd0));

endmodule
