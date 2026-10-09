// The paired refutation for unpacked_array_port.sv's output half: a
// memory driven by a submodule holds what that submodule wrote, so a
// different expected value must be refuted. Without this, an array
// the child never actually reached would satisfy the holds case too.
module unpacked_array_port_fails (
    input logic clk,
    input logic rst,
    input logic [1:0] i,
    input logic [7:0] d
);

  logic [7:0] driven[0:3];

  array_writer2 u_write (
      .clk(clk),
      .rst(rst),
      .i  (i),
      .d  (d),
      .m  (driven)
  );

  assert property (@(posedge clk) (!rst && !$past(rst) && $past(i) == i)
                   |-> driven[i] == $past(d) + 8'd1);

endmodule

module array_writer2 (
    input logic clk,
    input logic rst,
    input logic [1:0] i,
    input logic [7:0] d,
    output logic [7:0] m[0:3]
);
  always_ff @(posedge clk) begin
    if (rst) m <= '{default: 0};
    else m[i] <= d;
  end
endmodule
