// Positional parameter override (`#(3)`) of parameter_counter
// (duplicated here from parameter.sv so this test is self-contained,
// since check_bmc() loads a single file with no filelist), instead of
// the named-override style (`#(.WIDTH(3))`) covered by
// named_param_override.sv.  With WIDTH overridden to 3, MAX becomes 7
// and the child's own `assert property (count != MAX)` falsifies at
// cycle 8; if the override were ignored (falling back to
// parameter_counter's own default WIDTH=4, MAX=15), count could never
// reach 15 within 8 cycles and BMC would return UNKNOWN instead.
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

module positional_param_override (input logic clk, input logic rst);

  parameter_counter #(3) pc (
      .clk(clk),
      .rst(rst)
  );

endmodule
