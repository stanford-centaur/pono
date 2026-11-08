// `get_coverage()` returns the coverage number as a real (LRM 19.8,
// Table 19-5), which is the one way a covergroup stops being a pure
// observer: branch on it and the coverage state a simulator
// accumulates decides what the design does.
module covergroup_coverage_read (
    input logic clk,
    input logic [3:0] a
);

  logic [3:0] x;

  initial x = 4'd0;

  covergroup cg @(posedge clk);
    coverpoint a;
  endgroup

  cg cg_inst = new();

  always_ff @(posedge clk) begin
    if (cg_inst.get_coverage() > 50.0) x <= 4'd5;
  end

  assert property (@(posedge clk) x == 4'd0);

endmodule
