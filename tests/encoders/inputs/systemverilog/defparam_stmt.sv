// Out-of-scope: legacy `defparam` parameter override, as opposed to
// the `#(...)` styles covered in test_systemverilog_hierarchy.cpp.
module leaf #(
    parameter WIDTH = 4
) (
    input logic clk,
    input logic rst
);
  logic [WIDTH-1:0] count;
  always_ff @(posedge clk) begin
    if (rst) count <= '0;
    else count <= count + 1;
  end
endmodule

module defparam_stmt (input logic clk, input logic rst);
  leaf l (
      clk,
      rst
  );
  defparam l.WIDTH = 8;
  // With the override correctly applied, `l.count` (8 bits) can't wrap
  // back to 0 within 20 post-reset cycles; if the override is dropped
  // (this gap), `l.count` stays 4 bits wide and wraps at cycle 16.
  assert property (@(posedge clk) l.count != 0);
endmodule
