// A pattern case that is not exhaustive: no `default`, and the one
// item's pattern can fail to match. `y` is then unassigned on some
// paths through a combinational block, so it holds its value --
// the latch real synthesis infers, which the encoder models as a
// next-state hold. That is why the selector that decides this
// cycle's value is last cycle's.
module pattern_case_incomplete (input logic clk, input logic [3:0] x);
  logic [3:0] y;
  logic started;
  initial started = 1'b0;
  always_ff @(posedge clk) started <= 1'b1;
  always_comb begin
    case (x) matches
      4'd1 : y = 4'd10;
    endcase
  end
  assert property (@(posedge clk)
                   (started && $past(x) != 4'd1) |-> y == $past(y));
endmodule
