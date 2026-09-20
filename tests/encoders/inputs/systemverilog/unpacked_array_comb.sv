// An unpacked-array element write from a combinational block.  Only
// clocked (always_ff) array writes are modelled; leaving this one out
// would silently leave `lut` unconstrained, and an unconstrained array
// reads as anything -- so it is rejected instead.
module unpacked_array_comb (
    input logic clk,
    input logic [1:0] idx,
    input logic [7:0] d,
    output logic [7:0] o
);

  logic [7:0] lut[0:3];

  always_comb lut[idx] = d;
  assign o = lut[0];

  assert property (@(posedge clk) o == 8'd0);

endmodule
