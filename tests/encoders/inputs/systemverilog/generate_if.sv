// Generate-`if`: the (default) parameter selects between two
// differently-stepping always_ff blocks at elaboration time.  With
// the default USE_FAST=1, the *fast* branch (+2/cycle) must be the
// one actually elaborated; count(k) = 2*(k-1), so count == 6 first at
// cycle 4.  If generate-if incorrectly fell through to the slow
// (+1/cycle) branch, count wouldn't reach 6 until cycle 7, so this
// bound cleanly distinguishes the two.
module generate_if_counter #(
    parameter USE_FAST = 1
) (
    input logic clk,
    input logic rst
);

  logic [3:0] count;

  if (USE_FAST) begin : g
    always_ff @(posedge clk) begin
      if (rst) count <= 4'd0;
      else count <= count + 4'd2;
    end
  end else begin : g
    always_ff @(posedge clk) begin
      if (rst) count <= 4'd0;
      else count <= count + 4'd1;
    end
  end

  assert property (@(posedge clk) count != 4'd6);

endmodule
