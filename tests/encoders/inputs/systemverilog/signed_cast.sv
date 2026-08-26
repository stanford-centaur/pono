// `$signed()`/`$unsigned()` reinterpret an unsigned expression's bits
// as signed (or vice versa) without changing its width; a subsequent
// comparison using that cast must use the *cast's* signedness, not
// the original operand's. `a` is an unsigned 4-bit input: as raw bits
// 4'b1111 == 15 (unsigned), but reinterpreted via $signed it is -1,
// which must compare as less than 0.
module signed_cast (input logic clk, input logic [3:0] a);

  assert property (@(posedge clk) (a != 4'b1111) || ($signed(a) < 0));

endmodule
