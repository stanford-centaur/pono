// A `ref` argument aliases the caller's variable for the whole call
// rather than being copied back at the end, so the body's reads and
// writes interleave with the caller's -- not what inlining a copy-in
// copy-out models.
module task_call_ref_arg (
    input logic clk,
    input logic [3:0] a
);

  logic [3:0] b;

  task automatic viaref(input logic [3:0] x, ref logic [3:0] y);
    y = x;
  endtask

  always_ff @(posedge clk) viaref(a, b);

  assert property (@(posedge clk) b == b);

endmodule
