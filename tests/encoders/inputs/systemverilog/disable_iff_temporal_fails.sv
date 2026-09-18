// disable_iff_temporal.sv without the `disable iff (rst)`: the reset
// cycle is no longer exempt, so `nexttime started` is evaluated there
// and fails, because `started` is still low one cycle later.  Confirms
// the holds variant passes because of the exemption and not because
// the property was trivially true.
module disable_iff_temporal_fails (
    input logic clk,
    input logic rst
);

  logic started;

  always_ff @(posedge clk) begin
    if (rst) started <= 1'b0;
    else started <= 1'b1;
  end

  assert property (@(posedge clk) nexttime started);

endmodule
