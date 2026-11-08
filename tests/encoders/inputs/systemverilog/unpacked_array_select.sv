// Selecting between two whole unpacked arrays with a conditional
// expression. Both branches are array terms, so this is an Ite on the
// arrays themselves; asking either branch for a bit width would abort
// inside the solver rather than raise a PonoException.
module unpacked_array_select (
    input logic clk,
    input logic rst,
    input logic c
);

  logic [7:0] a[0:3];
  logic [7:0] b[0:3];
  logic [7:0] d[0:3];

  always_ff @(posedge clk) begin
    foreach (a[i]) a[i] <= 8'd1;
    foreach (b[i]) b[i] <= 8'd2;
    d <= c ? a : b;
  end

  assert property (@(posedge clk) (!rst && !$past(rst) && !$past(rst, 2))
                   |-> d[0] == ($past(c) ? 8'd1 : 8'd2));

endmodule
