// An unpacked array is passed whole: the child's port variable and
// the parent's array are the same term. Connecting a *slice*
// therefore cannot work, and must be refused rather than binding
// the parent's array to the port's sort -- which left the parent
// indexing its own four-element array with the port's one-bit
// index, and reached the solver as a sort mismatch rather than an
// error naming the port.
module unpacked_array_port_slice_child (
    output logic [7:0] o[0:1]
);
  assign o[0] = 8'd1;
  assign o[1] = 8'd2;
endmodule

module unpacked_array_port_slice (
    input logic clk,
    input logic rst
);

  logic [7:0] p[0:3];
  unpacked_array_port_slice_child u (.o(p[0:1]));

  assert property (@(posedge clk) p[0] == 8'd1);

endmodule
