// Indices outside an unpacked array's declared range. The LRM
// ignores such a write and reads X, and neither may disturb a real
// cell -- the normalized index is truncated to the array's address
// width, so without a range check `m[23]` of an `m[3:18]` folds onto
// `m[7]`.
//
// `keep` pins that: it is written through the in-range index 7 in
// the same cycle as the out-of-range 23, and must survive. `phantom`
// pins the other half: index 6 of a five-element array is read back
// after being written, which must *not* hold, so this fixture is
// paired with unpacked_array_out_of_range_fails.sv rather than
// asserting it here.
module unpacked_array_out_of_range (
    input logic clk,
    input logic rst
);

  logic [7:0] m[3:18];
  logic [7:0] keep;

  always_ff @(posedge clk) begin
    m[7]  <= 8'd11;
    m[23] <= 8'd99;
    keep  <= m[7];
  end

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> keep == 8'd11);

endmodule
