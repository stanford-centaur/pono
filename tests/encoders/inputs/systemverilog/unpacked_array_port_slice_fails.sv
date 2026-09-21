// Elements outside the slice are not driven, and the slice's offset
// is real: `p[0]` is neither driven nor the child's first element.
// Must be refuted.
module unpacked_array_port_slice_fails_child (
    output logic [7:0] o[0:1]
);
  assign o[0] = 8'd1;
  assign o[1] = 8'd2;
endmodule

module unpacked_array_port_slice_fails (
    input logic clk,
    input logic rst
);

  logic [7:0] p[0:3];
  unpacked_array_port_slice_fails_child u (.o(p[1:2]));

  assert property (@(posedge clk) p[0] == 8'd1);

endmodule
