// `priority if` and `unique case` as semantic modifiers (not just a
// parse-through of plain if/case): two independent decision
// structures, each driven by a free selector, combined into one
// property so a single falsifying trace has to pick both selectors
// correctly in the same cycle.  Paired with
// priority_if_unique_case_holds.sv.
module priority_if_unique_case (
    input logic clk,
    input logic rst,
    input logic [1:0] s,
    input logic sel
);

  logic [3:0] x, y;

  always_ff @(posedge clk) begin
    if (rst) begin
      x <= 4'd0;
      y <= 4'd0;
    end else begin
      priority if (sel) x <= 4'd5;
      else x <= 4'd3;

      unique case (s)
        2'd0: y <= 4'd1;
        2'd1: y <= 4'd2;
        2'd2: y <= 4'd4;
        2'd3: y <= 4'd8;
      endcase
    end
  end

  assert property (@(posedge clk) x != 4'd5 || y != 4'd8);

endmodule
