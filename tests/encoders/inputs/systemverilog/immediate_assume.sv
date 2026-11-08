// Procedural immediate `assume` (distinct from concurrent `assume
// property (...)`, covered by assume_property.sv): ImmediateAssertion
// handles Assume/Restrict the same way it handles Assert. `data` is
// otherwise a completely free input, so if the assumption weren't
// honored, BMC could pick data == 7 freely.
module immediate_assume (input logic clk, input logic rst, input logic [3:0] data);

  always_comb begin
    assume (data == 4'd0);
  end

  assert property (@(posedge clk) data != 4'd7);

endmodule
