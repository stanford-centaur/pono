// A count above one has to nest: reaching the second occurrence is
// reaching the first and then another after it. Exactly one grant
// ever arrives, so a second one never does and this is refuted --
// which a count that collapsed to one occurrence would not be.
module goto_repetition_count (
    input logic clk
);

  logic seen, gnt;

  initial seen = 1'b0;

  always_ff @(posedge clk) begin
    gnt  <= !seen;
    seen <= 1'b1;
  end

  assert property (@(posedge clk) gnt[->2]);

endmodule
