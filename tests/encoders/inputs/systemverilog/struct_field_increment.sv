// `++` on a packed-struct field, not just a plain scalar.
// p.cnt(k) = k-1 for k >= 1, so p.cnt first equals 3 at cycle 4.
typedef struct packed {
  logic [3:0] cnt;
  logic       flag;
} pair_t;

module struct_field_increment (input logic clk, input logic rst);

  pair_t p;

  always_ff @(posedge clk) begin
    if (rst) p.cnt <= 4'd0;
    else p.cnt++;
  end

  assert property (@(posedge clk) p.cnt != 4'd3);

endmodule
