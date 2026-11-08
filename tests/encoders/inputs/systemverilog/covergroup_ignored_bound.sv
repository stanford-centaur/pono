// The same design and covergroup with a bound the counter passes,
// so the covergroup is shown not to suppress a counterexample
// either.
module covergroup_ignored_bound (
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

  assert property (@(posedge clk) cnt <= 4'd5);

endmodule
