// The same declarations with a bound the counter passes, so they
// are shown not to suppress a counterexample either.
module unrepresentable_decls_bound (
    input logic clk,
    input logic rst
);

  // None of these is referenced anywhere below.
  event   done;
  real    threshold;
  string  label;
  chandle ctx;

  logic [3:0] cnt;

  initial cnt = 4'd0;

  always_ff @(posedge clk) begin
    if (rst) cnt <= 4'd0;
    else if (cnt != 4'd9) cnt <= cnt + 4'd1;
  end

  assert property (@(posedge clk) cnt <= 4'd5);

endmodule
