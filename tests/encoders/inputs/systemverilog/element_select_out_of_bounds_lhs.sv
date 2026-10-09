// A constant element-select lvalue whose index is out of range for
// its base (`flag[10]` into a 4-bit `flag`). The LRM makes such a
// write a no-op, and that now falls out rather than being special-
// cased: the index names no bit range, so the lvalue resolver
// declines and the write takes the dynamic-position splice, which
// masks by the shifted width -- and a position past the end shifts
// that mask away entirely, leaving the value untouched.
//
// `flag` therefore has no effective driver and holds whatever it
// started with. `started` gates the first cycle, where `$past` has
// no history to compare against.
module element_select_out_of_bounds_lhs (
    input logic clk,
    input logic a
);

  logic [3:0] flag;
  logic started;

  initial started = 1'b0;

  always_ff @(posedge clk) begin
    flag[10] <= a;
    started  <= 1'b1;
  end

  assert property (@(posedge clk) started |-> flag == $past(flag));

endmodule
