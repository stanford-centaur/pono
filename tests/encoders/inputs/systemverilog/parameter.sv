module parameter_counter #(parameter int WIDTH = 4)
                          (input logic clk, input logic rst);

  localparam logic [WIDTH-1:0] MAX = (1 << WIDTH) - 1;

  logic [WIDTH-1:0] count;

  always @(posedge clk) begin
    if (rst)               count <= '0;
    else if (count != MAX) count <= count + 1;
    else                   count <= '0;
  end

  assert property (@(posedge clk) count != MAX);

endmodule
