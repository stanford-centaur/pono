// An in-property `case` (`AssertionExprKind::Case`) as a top-level
// property connective, analogous to property_conditional.sv's `if`.
module property_case (
    input logic clk,
    input logic [1:0] sel,
    input logic a,
    input logic b
);
  assert property (@(posedge clk) case (sel)
    0 : a;
    1 : b;
    default : 1'b1;
  endcase);
endmodule
