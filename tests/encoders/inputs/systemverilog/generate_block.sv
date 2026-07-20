// Generate-for block creating four independent counters.  Each
// iteration's body lives in its own `ctr[i]` scope; the genvar
// `i` is exposed by slang as a per-iteration ParameterSymbol.
// The assertion outside the loop reads one of the per-iteration
// counters via a hierarchical reference (HierarchicalValue), so
// this test also covers genvar -> parameter, generate-block
// member walking, and cross-scope reads.
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
