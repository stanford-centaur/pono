// Temporaries declared inside a procedural block. `t` stages a
// computation, `pick` is overwritten on one path, and `k` is a
// constant the unrolling machinery still folds. None of them is a
// register: each is bound to the value it carries at that point in
// the block, so reads of it mean that value.
//
// An uninitialized 4-state local defaults to all-X, which has no
// counterpart in this 2-valued model -- the declaration must bind
// nothing rather than invent a number, which is what the assertions
// below would catch.
module procedural_temp (
    input logic clk,
    input logic rst,
    input logic c,
    input logic [7:0] a,
    input logic [7:0] b
);

  logic [8:0] sum;
  logic [7:0] chosen;
  logic [7:0] scaled;

  always_ff @(posedge clk) begin
    logic [8:0] t;
    logic [7:0] pick;
    int k = 3;
    t = a + b;
    pick = a;
    if (c) pick = b;
    sum <= t;
    chosen <= pick;
    scaled <= a + 8'(k);
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (sum == $past(a) + $past(b)
                        && chosen == ($past(c) ? $past(b) : $past(a))
                        && scaled == $past(a) + 8'd3));

endmodule
