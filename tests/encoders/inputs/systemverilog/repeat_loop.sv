// `repeat` version of while_loop.sv: a fixed trip count with a
// manually-managed index rather than a loop-carried condition.
module repeat_loop (input logic clk, input logic rst, input logic [3:0] din);

  logic any_set, reg_any_set;

  always_comb begin
    int i;
    any_set = 1'b0;
    i = 0;
    repeat (4) begin
      any_set |= din[i];
      i = i + 1;
    end
  end

  always_ff @(posedge clk) begin
    if (rst) reg_any_set <= 1'b0;
    else reg_any_set <= any_set;
  end

  assert property (@(posedge clk) reg_any_set != 1'b1);

endmodule
