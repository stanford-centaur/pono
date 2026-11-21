// case_inside_range.sv's range widened by one at the top end. Must
// be refuted, or the range item would not be bounding anything.
module case_inside_range_fails (
    input logic clk,
    input logic rst,
    input logic [3:0] s
);

  logic yr;

  always_ff @(posedge clk) begin
    case (s) inside
      [4'd1:4'd3]: yr <= 1'b1;
      default: yr <= 1'b0;
    endcase
  end

  assert property (@(posedge clk) (!rst && !$past(rst))
                   |-> (yr == $past(s >= 4'd1 && s <= 4'd4)));

endmodule
