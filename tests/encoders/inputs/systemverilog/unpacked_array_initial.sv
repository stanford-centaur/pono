// Unpacked-array elements written from an `initial` block, which
// constrains the initial state rather than a next one. The array
// holds thereafter, so the values have to survive to be read.
module unpacked_array_initial (
    input logic clk
);

  logic [7:0] m[0:3];

  initial begin
    m[1] = 8'd7;
    m[2] = 8'd8;
  end

  always_ff @(posedge clk) m <= m;

  assert property (@(posedge clk) (m[1] == 8'd7 && m[2] == 8'd8));

endmodule
