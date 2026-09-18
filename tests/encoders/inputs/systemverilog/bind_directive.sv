// Out-of-scope: the `bind` directive, instantiating a checker-like
// module into another module's scope from outside its source.
module leaf2 (input logic clk, input logic rst);
  logic [3:0] count;
  always_ff @(posedge clk) begin
    if (rst) count <= '0;
    else count <= count + 1;
  end
endmodule

module checker_mod (input logic clk, input logic [3:0] count);
  // References leaf2's own internal `count` directly -- bind's whole
  // point is attaching an assertion that sees the target's internals
  // without modifying its source. A small bound is deliberately used
  // (count can be at most 3 after 3 post-reset cycles) so a
  // disconnected/free `count` (bind's port binding silently failing)
  // would trivially falsify this immediately, distinguishing that
  // failure mode from "bind isn't processed at all".
  assert property (@(posedge clk) count != 4'd15);
endmodule

module bind_directive (input logic clk, input logic rst);
  leaf2 l2 (
      clk,
      rst
  );
  bind leaf2 checker_mod cm (.clk(clk), .count(count));
endmodule
