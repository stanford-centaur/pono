// Generate-`case`: the parameter selects among three differently-
// stepping always_ff blocks at elaboration time. With the default
// MODE=1 (+2/cycle), count == 6 first at cycle 4; MODE=0 (+1/cycle)
// wouldn't reach 6 until cycle 7, and the default branch (+3/cycle)
// would overshoot 6 without hitting it exactly -- so the bound
// identifies which branch actually elaborated.
module generate_case_counter #(
    parameter int MODE = 1
) (
    input logic clk,
    input logic rst
);

  logic [3:0] count;

  case (MODE)
    0: begin : g
      always_ff @(posedge clk) begin
        if (rst) count <= 4'd0;
        else count <= count + 4'd1;
      end
    end
    1: begin : g
      always_ff @(posedge clk) begin
        if (rst) count <= 4'd0;
        else count <= count + 4'd2;
      end
    end
    default: begin : g
      always_ff @(posedge clk) begin
        if (rst) count <= 4'd0;
        else count <= count + 4'd3;
      end
    end
  endcase

  assert property (@(posedge clk) count != 4'd6);

endmodule
