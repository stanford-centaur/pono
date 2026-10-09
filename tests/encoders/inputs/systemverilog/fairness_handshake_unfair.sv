// fairness_handshake.sv with nothing promising the acknowledgement:
// `ack` is free and may never arrive, leaving `busy` set forever.
module fairness_handshake_unfair (
    input logic clk,
    input logic req,
    input logic ack
);

  logic busy;
  initial busy = 1'b0;
  always_ff @(posedge clk) begin
    if (ack) busy <= 1'b0;
    else if (req) busy <= 1'b1;
  end

  assert property (@(posedge clk) s_eventually !busy);

endmodule
