// A streaming concatenation as an output-port connection. `>>`
// re-orders nothing, so the port's bits land on the stream's
// expressions exactly as a plain concatenation connection would --
// one output-alias segment per expression. (`<<` moves bits across
// those boundaries, so it needs a segment per block instead -- see
// streaming_port_lshift.sv.)
module streaming_concat_port_child (
    output logic [7:0] o
);

  assign o = 8'hA5;

endmodule

module streaming_concat_port (
    input logic clk
);

  logic [3:0] hi, lo;

  streaming_concat_port_child u (.o({>>{hi, lo}}));

  assert property (@(posedge clk) hi == 4'hA && lo == 4'h5);

endmodule
