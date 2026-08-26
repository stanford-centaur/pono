// `final` blocks run once at the end of simulation for cleanup/
// reporting ($display, checks) -- they have no synthesis meaning and
// no analog in an infinite/bounded-trace model-checking model (there
// is no "end of simulation" here), so this should be intentionally
// ignored, the same as $display and other simulation-only constructs
// elsewhere in this encoder, not something that needs modeling.
module final_block (input logic clk, input logic rst, input logic [3:0] data);

  logic [3:0] data_reg;

  always_ff @(posedge clk) begin
    if (rst) data_reg <= 4'd0;
    else data_reg <= data;
  end

  final begin
    $display("final data_reg = %d", data_reg);
  end

  assert property (@(posedge clk) data_reg == data_reg);

endmodule
