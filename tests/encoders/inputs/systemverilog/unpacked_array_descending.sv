// A constant assignment pattern over a *descending* unpacked range.
// A pattern lists its values left to right, which is the order of
// the declared indices -- so `'{10, 20, 30, 40}` into `m[3:0]` puts
// 10 at m[3] and 40 at m[0]. Position therefore counts the opposite
// way from the normalized index the array sort uses, and reading it
// straight across fills the array backwards.
//
// Element-by-element writes and `foreach` were never affected, only
// the pattern, so both are asserted here too: the reversal would
// otherwise be invisible whenever two of them cancel out.
module unpacked_array_descending (
    input logic clk,
    input logic rst
);

  logic [7:0] m[3:0];
  logic [7:0] n[3:0];
  logic [7:0] p[0:3];

  always_ff @(posedge clk) begin
    m <= '{8'd10, 8'd20, 8'd30, 8'd40};
    foreach (n[i]) n[i] <= 8'(i);
    p <= '{8'd10, 8'd20, 8'd30, 8'd40};
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (m[3] == 8'd10 && m[2] == 8'd20
                        && m[1] == 8'd30 && m[0] == 8'd40
                        && n[3] == 8'd3 && n[0] == 8'd0
                        && p[0] == 8'd10 && p[3] == 8'd40));

endmodule
