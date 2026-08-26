// Generate-for block creating four independent counters, each in its
// own `ctr[i]` scope. The assertion outside the loop reads one
// iteration's counter via a hierarchical reference (`ctr[2].count`),
// covering generate-block member walking and cross-scope reads.
module generate_block_counters (input logic clk, input logic rst);

  for (genvar i = 0; i < 4; ++i) begin : ctr
    logic [3:0] count;
    always_ff @(posedge clk) begin
      if (rst) count <= '0;
      else     count <= count + 1;
    end
  end

  // ctr[2].count reaches 5 at cycle 6 (one reset cycle + five
  // increments from 0).
  assert property (@(posedge clk) ctr[2].count != 5);

endmodule
