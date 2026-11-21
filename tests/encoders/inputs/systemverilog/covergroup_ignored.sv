// Functional coverage observes the design and never drives it, so
// a covergroup can be walked past: the counter's bound is proved
// here with one in scope.
module covergroup_ignored (
    input logic clk,
    input logic rst
);

  logic [3:0] cnt;

  initial cnt = 4'd0;

  always_ff @(posedge clk) begin
    if (rst) cnt <= 4'd0;
    else if (cnt != 4'd9) cnt <= cnt + 4'd1;
  end

  covergroup cg @(posedge clk);
    coverpoint cnt;
  endgroup

  cg cg_inst = new();

  assert property (@(posedge clk) cnt <= 4'd9);

endmodule
