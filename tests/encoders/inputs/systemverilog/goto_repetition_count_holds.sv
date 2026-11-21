// The other side of goto_repetition_count.sv: grants arrive forever,
// so any number of occurrences is eventually reached.
module goto_repetition_count_holds (
    input logic clk
);

  logic gnt;

  always_ff @(posedge clk) gnt <= !gnt;

  assert property (@(posedge clk) gnt[->2]);

endmodule
