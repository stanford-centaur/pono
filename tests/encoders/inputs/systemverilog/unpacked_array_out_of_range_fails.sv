// Index 6 does not exist in a five-element array, so the write is
// ignored and the read is X: this read-back must be refuted. It was
// proved before the range check, certifying a cell the array does
// not have.
module unpacked_array_out_of_range_fails (
    input logic clk
);

  logic [7:0] m[0:4];
  logic [7:0] rd;

  always_ff @(posedge clk) begin
    m[3'd6] <= 8'd42;
    rd <= m[3'd6];
  end

  assert property (@(posedge clk) rd == 8'd42);

endmodule
