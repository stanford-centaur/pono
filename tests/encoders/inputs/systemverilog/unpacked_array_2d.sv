// A multi-dimensional unpacked array.  Only a single dimension with an
// integral element type is modelled, so this is rejected rather than
// silently losing a dimension.
module unpacked_array_2d (
    input logic clk,
    input logic rst,
    input logic [7:0] wdata
);

  logic [7:0] mem[0:3][0:3];

  always_ff @(posedge clk) begin
    if (rst) mem[0][0] <= 8'd0;
    else mem[0][0] <= wdata;
  end

  assert property (@(posedge clk) mem[0][0] == $past(wdata));

endmodule
