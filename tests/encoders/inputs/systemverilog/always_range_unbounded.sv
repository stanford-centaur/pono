// `always [m:$]` -- the one windowed form whose upper bound may be
// unbounded, encoded as m forward shifts around the ordinary G tester.
// `started` is cleared while rst is high and latches high afterwards,
// so from cycle 2 onwards it is high forever.
module always_range_unbounded (
    input logic clk,
    input logic rst
);

  logic started;

  always_ff @(posedge clk) begin
    if (rst) started <= 1'b0;
    else started <= 1'b1;
  end

  assert property (@(posedge clk) always [2:$] started);

endmodule
