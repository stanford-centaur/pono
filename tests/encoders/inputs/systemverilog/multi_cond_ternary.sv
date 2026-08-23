// A ternary condition can be a `&&&`-joined list of more than one
// expression (LRM 11.4.11), legal outside case/if context too --
// `a &&& b` means both must hold. expr_to_term() previously only ever
// read conditions[0], silently ignoring `&&& b`: with a true but b
// false, the buggy encoding picks the "true" branch anyway.
module multi_cond_ternary (input logic clk,
                           input logic a,
                           input logic b);

  logic [3:0] result;

  assign result = a &&& b ? 4'd1 : 4'd0;

  assert property (@(posedge clk) !(a && !b) || (result == 4'd0));

endmodule
