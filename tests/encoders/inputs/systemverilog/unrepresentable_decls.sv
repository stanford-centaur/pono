// Variables whose types have no bit-vector representation,
// declared and never used. A declaration computes nothing, so
// these leave the rest of the module exactly as it would be
// without them -- the counter's bound is still proved.
module unrepresentable_decls (
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

  assert property (@(posedge clk) cnt <= 4'd9);

endmodule
