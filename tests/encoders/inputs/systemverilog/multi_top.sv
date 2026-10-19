// Two modules, neither instantiating the other, so slang makes both
// top-level -- and sorts the tops alphabetically by name, so
// `aaa_helper` comes first regardless of source order or which file
// was named on the command line.
//
// `aaa_helper`'s assertion is trivially true; `dut`'s is genuinely
// false.  Taking the first top without saying so therefore reports a
// clean pass for a design containing a violated assertion, which is
// why an ambiguous top has to be rejected instead.
module aaa_helper (
    input logic clk,
    input logic rst
);

  logic [1:0] c;

  always_ff @(posedge clk) begin
    if (rst) c <= 2'd0;
    else c <= c + 2'd1;
  end

  // A 2-bit counter is always <= 3.
  assert property (@(posedge clk) c <= 2'd3);

endmodule

module dut (
    input logic clk,
    input logic rst
);

  logic [1:0] c;

  always_ff @(posedge clk) begin
    if (rst) c <= 2'd0;
    else c <= c + 2'd1;
  end

  // False: the counter does reach 3.
  assert property (@(posedge clk) c != 2'd3);

endmodule
