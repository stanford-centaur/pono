// `break` inside the already-supported `for` loop.  `din` is a fixed
// constant (4'b1000, only bit 3 set) rather than a free input so the
// outcome unambiguously distinguishes "break honored" from "break
// ignored": correct semantics stop the loop before i reaches 2, so
// bit 3 is never visited and any_set stays 0 forever (property
// holds); if break were silently ignored, the loop would reach bit 3
// and any_set would go high (falsifying at cycle 2).
module break_in_for (input logic clk, input logic rst);

  logic [3:0] din;
  logic       any_set, reg_any_set;

  assign din = 4'b1000;

  always_comb begin
    any_set = 1'b0;
    for (int i = 0; i < 4; ++i) begin
      if (i == 2) break;
      any_set |= din[i];
    end
  end

  always_ff @(posedge clk) begin
    if (rst) reg_any_set <= 1'b0;
    else reg_any_set <= any_set;
  end

  assert property (@(posedge clk) reg_any_set != 1'b1);

endmodule
