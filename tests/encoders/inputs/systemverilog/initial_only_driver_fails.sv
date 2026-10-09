// The other half of initial_only_driver.sv: an element the initial
// block never wrote is not zero, it is unknown. Claiming a value
// for it must be refuted -- holding the *whole* array at its
// initial value would let this pass.
module initial_only_driver_fails (
    input logic clk,
    input logic rst
);

  logic [7:0] m[0:3];

  initial m[1] = 8'd7;

  assert property (@(posedge clk) m[2] == 8'd0);

endmodule
