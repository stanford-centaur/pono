// Case equality (`===`) and wildcard equality (`==?`) registered
// together: `a === 4'b1010` exactly implies `a ==? 4'b10??` (the
// wildcard pattern's don't-care bits already match 1010's low two
// bits), so the implication genuinely holds for every a.
module equality_variants (input logic clk, input logic rst, input logic [3:0] a);

  logic case_eq, wild_eq;

  always_ff @(posedge clk) begin
    if (rst) begin
      case_eq <= 1'b0;
      wild_eq <= 1'b0;
    end else begin
      case_eq <= (a === 4'b1010);
      wild_eq <= (a ==? 4'b10??);
    end
  end

  assert property (@(posedge clk) case_eq |-> wild_eq);

endmodule
