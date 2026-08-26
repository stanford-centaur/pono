// `++` on a packed-array element (constant index), not just a plain
// scalar.
// mem[1](k) = k-1 for k >= 1, so mem[1] first equals 3 at cycle 4.
module element_increment (input logic clk, input logic rst);

  logic [3:0][3:0] mem;

  always_ff @(posedge clk) begin
    if (rst) mem <= '0;
    else mem[1]++;
  end

  assert property (@(posedge clk) mem[1] != 4'd3);

endmodule
