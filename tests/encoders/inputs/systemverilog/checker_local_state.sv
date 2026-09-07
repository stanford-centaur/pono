// A checker with its own local sequential state (`count`, driven by
// an always_ff inside the checker body) rather than only referencing
// its formal ports -- exercises SymbolTable::pre_scan_state_vars()'s
// recursion into checker instances (to classify `count` as a state
// var before anything is declared) and InstanceEncoder::
// process_checker_instance()'s own declare pass for it.
checker my_checker (input logic clk, input logic rst, input logic a);
  logic [3:0] count;
  always_ff @(posedge clk) begin
    if (rst) count <= 4'd0;
    else if (a) count <= count + 4'd1;
  end
  assert property (@(posedge clk) count < 4'd3);
endchecker

module checker_local_state (
    input logic clk,
    input logic rst,
    input logic a
);
  my_checker chk (clk, rst, a);
endmodule
