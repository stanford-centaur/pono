// An assignment pattern whose values are not elaboration-time
// constants. The constant case folds to a single constant array;
// this one has a term per element, so it builds as stores over a
// seed -- seeded with the first value, which keeps a uniform
// pattern (the common case) a constant array with no stores at all.
//
// `m` is ascending and `d` descending, and the pattern is
// deliberately non-uniform, since position-to-index ordering is
// invisible when every element is the same. A pattern lists its
// values in declared-index order, so `'{v, w, w, w}` puts v at m[0]
// and at d[3].
module unpacked_pattern_nonconstant (
    input logic clk,
    input logic rst,
    input logic [7:0] v,
    input logic [7:0] w
);

  logic [7:0] m[0:3];
  logic [7:0] d[3:0];
  logic [7:0] x[0:3];

  always_ff @(posedge clk) begin
    m <= '{v, w, w, w};
    d <= '{v, w, w, w};
    x <= '{8'd7, v, v, v};
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (m[0] == $past(v) && m[1] == $past(w)
                        && m[3] == $past(w)
                        && d[3] == $past(v) && d[0] == $past(w)
                        && x[0] == 8'd7 && x[1] == $past(v)));

endmodule
