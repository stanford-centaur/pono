// The paired refutation for unbounded_delay.sv: the grant never
// comes, so a request that is never answered is a fair lasso
// violating the eventuality. Without this, an F whose justice
// condition is never emitted -- so never has to be discharged --
// would satisfy the holds case just as well.
module unbounded_delay_fails (
    input logic clk,
    input logic req
);

  logic gnt;

  always_ff @(posedge clk) gnt <= 1'b0;

  assert property (@(posedge clk) req |-> ##[1:$] gnt);

endmodule
