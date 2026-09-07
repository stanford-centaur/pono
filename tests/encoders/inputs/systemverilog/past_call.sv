// `$past(expr)` returns expr's value from the previous cycle.
// The encoder implements this via a hidden 1-cycle latch state
// var.  The register `reg_data` already tracks the previous
// `data`, so the property `reg_data != $past(data)` is falsified
// once they agree (which happens immediately after reset releases).
module past_call (input logic clk, input logic rst,
                  input logic [3:0] data);
  logic [3:0] reg_data;
  always_ff @(posedge clk) reg_data <= rst ? '0 : data;
  assert property (@(posedge clk) reg_data != $past(data));
endmodule
