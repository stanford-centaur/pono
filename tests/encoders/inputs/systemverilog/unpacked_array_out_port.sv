// An unpacked array driven through a child's output port, with the
// child driving its own port array element by element -- which is
// the same continuous-assign-to-an-array-element path, reached from
// inside the child.
module unpacked_array_out_port_child (
    output logic [7:0] o[0:1]
);
  assign o[0] = 8'd1;
  assign o[1] = 8'd2;
endmodule

module unpacked_array_out_port (
    input logic clk,
    input logic rst
);

  logic [7:0] p[0:1];
  unpacked_array_out_port_child u (.o(p));

  assert property (@(posedge clk) p[0] == 8'd1 && p[1] == 8'd2);

endmodule
