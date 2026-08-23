// `case (x) matches ... endcase` (pattern-matching case,
// StatementKind::PatternCase) is a distinct statement kind from plain
// `case`/`casex`/`casez` (StatementKind::Case). process_statement()'s
// main switch has no explicit case for it, but its generic
// unhandled-statement-kind default already logs a warning and skips
// it cleanly (confirmed empirically) -- this locks that behavior in.
module pattern_case (input logic clk, input logic [3:0] x);

  logic [3:0] y;

  always_ff @(posedge clk) begin
    case (x) matches
      4'd1: y <= 4'd10;
      default: y <= 4'd0;
    endcase
  end

  assert property (@(posedge clk) 1'b1);

endmodule
