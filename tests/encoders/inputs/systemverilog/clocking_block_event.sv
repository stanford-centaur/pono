// A named clocking block used as an assertion's clocking event,
// `@(cb)`. slang binds that as a void-typed reference to the block
// rather than to a signal, so the clock is the event the block itself
// was declared with -- and every check the written-out form gets
// applies to that event too, since resolving it recurses.
//
// A `default clocking` block needs no such handling: slang
// substitutes its event during elaboration, so an assertion with no
// clocking event of its own arrives already carrying `@(posedge
// clk)`. This is the form that does not.
module clocking_block_event (
    input logic clk,
    input logic [3:0] d
);

  logic [3:0] q;
  logic started;

  clocking cb @(posedge clk);
  endclocking

  initial started = 1'b0;

  always_ff @(posedge clk) begin
    q <= d;
    started <= 1'b1;
  end

  assert property (@(cb) started |-> q == $past(d));

endmodule
