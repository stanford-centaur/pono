// Named parameter override (`#(.WIDTH(3))`) of parameter_counter
// (duplicated here from parameter.sv so this test is self-contained).
// Same distinguishing bound as positional_param_override.sv, checking
// the other override spelling.
module parameter_counter #(
    parameter int WIDTH = 4
) (
    input logic clk,
    input logic rst
);

  localparam logic [WIDTH-1:0] MAX = (1 << WIDTH) - 1;

  logic [WIDTH-1:0] count;

  always @(posedge clk) begin
    if (rst) count <= '0;
    else if (count != MAX) count <= count + 1;
    else count <= '0;
  end

  assert property (@(posedge clk) count != MAX);

endmodule

module named_param_override (input logic clk, input logic rst);

  parameter_counter #(
      .WIDTH(3)
  ) pc (
      .clk(clk),
      .rst(rst)
  );

endmodule
