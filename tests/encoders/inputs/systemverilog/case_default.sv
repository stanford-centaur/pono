// Regression for a case-statement bug where `default:` used to apply
// unconditionally alongside whichever other item already matched: if
// that bug were present, x would always become 9 (the default value)
// regardless of `s`, and `x != 1` would hold forever.  With the fix,
// `s == 0` lets x reach 1 at the first post-reset transition,
// falsifying the property at cycle 2.
module case_default (input logic clk, input logic rst, input logic [1:0] s);

  logic [3:0] x;

  always_ff @(posedge clk) begin
    if (rst) x <= 4'd0;
    else
      case (s)
        2'd0: x <= 4'd1;
        2'd1: x <= 4'd2;
        default: x <= 4'd9;
      endcase
  end

  assert property (@(posedge clk) x != 4'd1);

endmodule
