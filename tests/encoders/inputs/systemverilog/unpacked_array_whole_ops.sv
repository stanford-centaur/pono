// Copying one whole array into another and comparing two whole
// arrays: SMT arrays support both directly, so each is a single term
// rather than an element-by-element expansion.
//
// `a` is written only at reset and then holds, and `b` tracks `a`, so
// the two agree at every cycle after reset -- an invariant that is
// inductive, and so provable rather than merely unrefuted.
module unpacked_array_whole_ops (
    input logic clk,
    input logic rst
);

  logic [7:0] a[0:3];
  logic [7:0] b[0:3];

  always_ff @(posedge clk) begin
    if (rst) begin
      a <= '{default: 0};
      b <= '{default: 0};
    end else begin
      b <= a;
    end
  end

  assert property (@(posedge clk) !rst |-> a == b);

endmodule
