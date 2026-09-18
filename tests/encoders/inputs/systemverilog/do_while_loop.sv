// `do ... while` version of while_loop.sv.
module do_while_loop (input logic clk, input logic rst, input logic [3:0] din);

  logic any_set, reg_any_set;

  always_comb begin
    int i;
    any_set = 1'b0;
    i = 0;
    do begin
      any_set |= din[i];
      i = i + 1;
    end while (i < 4);
  end

  always_ff @(posedge clk) begin
    if (rst) reg_any_set <= 1'b0;
    else reg_any_set <= any_set;
  end

  assert property (@(posedge clk) reg_any_set != 1'b1);

endmodule
