// A legacy `defparam` parameter override, as opposed to the `#(...)`
// styles covered elsewhere in test_systemverilog_hierarchy.cpp.
//
// slang's elaborator applies the override before the encoder sees the
// design, the same way it splices in a `bind`, so nothing here needs
// its own handling. What this has to check is that the override really
// reached the model rather than the declared default: only an 8-bit
// counter ever equals 20, since a 4-bit one wraps at 16. A
// counterexample is therefore the override having been applied, and
// its depth says which cycle the counter got there.
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

module defparam_stmt (
    input logic clk,
    input logic rst
);

  leaf l (
      clk,
      rst
  );
  defparam l.WIDTH = 8;

  assert property (@(posedge clk) l.count != 8'd20);

endmodule
