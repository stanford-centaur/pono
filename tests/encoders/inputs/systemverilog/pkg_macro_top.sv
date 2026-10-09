// Compiler-macro (`define`/`ifdef`) preprocessing combined with a
// package `import` resolved across a --sv-filelist boundary:
// common_pkg.sv (providing `nibble_t`/`STEP`) is named only by
// pkg_filelist.f, exactly like the module-hierarchy case already
// covered by filelist_top.sv/filelist.f.
`define GUARD_WIDTH 4
`ifdef GUARD_WIDTH
`else
`define GUARD_WIDTH 1
`endif

import common_pkg::*;

module pkg_macro_top (input logic clk, input logic rst);

  nibble_t count;

  always_ff @(posedge clk) begin
    if (rst) count <= 4'd0;
    else count <= count + STEP;
  end

  // count(k) = STEP*(k-1) = 2*(k-1), so count == 6 first at cycle 4,
  // *if* the imported typedef/localparam resolved correctly.
  assert property (@(posedge clk) count != 4'd6);

endmodule
