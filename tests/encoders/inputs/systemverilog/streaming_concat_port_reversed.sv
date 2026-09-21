// A `<<` streaming concatenation as an output-port connection. The
// reversal moves bits across the boundaries between the stream's
// expressions, and an output alias is one contiguous port range per
// symbol, so there is no segment list that describes the result --
// hence a clean refusal rather than the `>>` case's split.
module streaming_concat_port_reversed_child (
    output logic [7:0] o
);

  assign o = 8'hA5;

endmodule

module streaming_concat_port_reversed (
    input logic clk
);

  logic [3:0] hi, lo;

  streaming_concat_port_reversed_child u (.o({<<4{hi, lo}}));

  assert property (@(posedge clk) 1'b1);

endmodule
