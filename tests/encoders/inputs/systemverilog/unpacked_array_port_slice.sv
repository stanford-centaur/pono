// An unpacked array connected to a *slice* of a parent array. The
// two cannot share one term, since they have different lengths, so
// the port gets an array of its own tied to the parent's element by
// element. The count is known at elaboration, so that is a fixed
// handful of equalities rather than anything dynamic.
//
// `p[1:2]` starts at 1 deliberately: an offset dropped somewhere
// would land the child's first element on p[0], which
// unpacked_array_port_slice_fails.sv rules out. `q` is a descending
// parent range, where normalized indices count up from the range's
// lower bound whichever way it is written.
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
  logic [7:0] q[3:0];

  unpacked_array_port_slice_child u (.o(p[1:2]));
  unpacked_array_port_slice_child v (.o(q[2:1]));

  assert property (@(posedge clk)
                   p[1] == 8'd1 && p[2] == 8'd2
                   && q[1] == 8'd1 && q[2] == 8'd2);

endmodule
