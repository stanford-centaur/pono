// Reading a temporary before anything writes it reads as X per the
// LRM, so it is rejected rather than given an arbitrary value.
module procedural_temp_read_first (
    input logic clk
);

  logic [7:0] y;

  always_ff @(posedge clk) begin
    logic [7:0] t;
    y <= t;
  end

  assert property (@(posedge clk) y == y);

endmodule
