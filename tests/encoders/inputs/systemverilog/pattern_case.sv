// `case (x) matches ... endcase` (pattern-matching case,
// StatementKind::PatternCase) is a distinct statement kind from plain
// `case`/`casex`/`casez` (StatementKind::Case), and its items are
// patterns rather than values to compare against.
//
// A constant pattern is the degenerate case, testing equality just as
// a plain case item would -- which is what makes it the right shape
// to pin the basic wiring with.
module pattern_case (
    input logic clk,
    input logic [3:0] x
);

  logic [3:0] y;

  always_ff @(posedge clk) begin
    case (x) matches
      4'd1: y <= 4'd10;
      default: y <= 4'd0;
    endcase
  end

  assert property (@(posedge clk) x == 4'd1 |=> y == 4'd10);

endmodule
