// The undriven half of unpacked_array_net.sv: `m[3]` has no driver,
// so it reads as an unknown rather than zero. Claiming a value for
// it must be refuted -- treating an undriven net as zero is the
// silent wrong answer this rules out.
module unpacked_array_net_fails (
    input logic clk,
    input logic rst,
    input logic [7:0] v
);

  wire [7:0] m[0:3];

  assign m[0] = v;

  assert property (@(posedge clk) m[3] == 8'd0);

endmodule
