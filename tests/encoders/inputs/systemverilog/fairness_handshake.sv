// What a fairness assumption is actually for. `busy` is set by a
// request and cleared only by an acknowledgement, so nothing in the
// design stops it staying set forever -- the environment has to
// promise the acknowledgement comes, which is what the assumption
// says.
//
// Paired with fairness_handshake_unfair.sv, the same design with no
// promise, which is refuted.
module fairness_handshake (
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

  assume property (@(posedge clk) s_eventually ack);
  assert property (@(posedge clk) s_eventually !busy);

endmodule
