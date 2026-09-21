// The paired refutation for goto_repetition.sv: the grant never
// comes, so the first occurrence never arrives.
module goto_repetition_fails (
    input logic clk,
    input logic req
);

  logic gnt;

  always_ff @(posedge clk) gnt <= 1'b0;

  assert property (@(posedge clk) req |-> gnt[->1]);

endmodule
