// Blocking `=` inside a clocked block. It infers a register exactly
// as `<=` does; the operator decides only when the new value becomes
// visible to later reads in the same block, which is immediately.
//
// The chain p/q/r therefore all land on b in one cycle, mix keeps the
// value a had *before* its blocking write only where `<=` is used,
// and the self-update reads the old n. Each conjunct is refuted by
// the opposite (non-blocking) reading, so this pins the semantics
// rather than merely exercising the path.
module blocking_in_clocked_block (
    input logic clk,
    input logic rst,
    input logic [7:0] b
);

  logic [7:0] p, q, r;
  logic [7:0] a, nb;
  logic [7:0] n;

  always_ff @(posedge clk) begin
    p = b;
    q = p;
    r = q;
  end

  always_ff @(posedge clk) begin
    a  = b;
    nb <= a;
  end

  always_ff @(posedge clk) begin
    if (rst) n = 8'd0;
    else n = n + 8'd1;
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (p == $past(b) && q == $past(b) && r == $past(b)
                        && nb == a && n == $past(n) + 8'd1));

endmodule
