// An unpacked array crossing a module boundary.  Supported unpacked
// arrays are registers internal to one module; a port would have to be
// spliced by the port-connection machinery, which is entirely
// bit-range arithmetic and has no array path.  Rejected rather than
// mis-encoded.
module unpacked_array_port (
    input logic clk,
    input logic [7:0] mem[0:3]
);

  assert property (@(posedge clk) mem[0] == 8'd0);

endmodule
