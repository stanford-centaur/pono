// A `$display` call used as a bare statement inside an `always_ff`
// block -- simulation-only output with no synthesis meaning, so
// process_statement()'s Call-expression handling recognizes it as a
// system call and skips it (logged), the same as the `$display`
// inside a `final` block in final_block.sv, rather than throwing the
// way a user-defined task call with real side effects does (see
// Gap_VoidTaskCall in test_systemverilog_statements.cpp).
module display_call_statement (
    input logic clk,
    input logic rst,
    input logic [3:0] a
);
  logic [3:0] b;
  always_ff @(posedge clk) begin
    if (rst) b <= 4'd0;
    else begin
      $display("a = %d", a);
      b <= a;
    end
  end
  assert property (@(posedge clk) rst || b == $past(a));
endmodule
