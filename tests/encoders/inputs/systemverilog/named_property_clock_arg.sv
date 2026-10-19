// A clocking-event argument naming a different clock than the rest
// of the design: the multiclock rejection must still fire.
module named_property_clock_arg (
    input logic clk,
    input logic clk2,
    input logic a,
    input logic b
);

  property p(ev, x);
    @(ev) x |-> ##1 x;
  endproperty
  assert property (@(posedge clk) a |-> ##1 b);
  assert property (p(posedge clk2, a));
endmodule
