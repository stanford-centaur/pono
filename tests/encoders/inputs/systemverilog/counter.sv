module counter (
    input logic clk,
    input logic rst
);

  logic [4:0] count;

  always @(posedge clk) begin
    if (rst) count <= 0;
    else count <= count + 1;
  end

  assert property (@(posedge clk) count != 4);

endmodule
