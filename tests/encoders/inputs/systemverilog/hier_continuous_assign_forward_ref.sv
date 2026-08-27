// A continuous assign whose target is a child instance's internal
// (non-port) signal, referenced via a hierarchical dot-path, where the
// assign appears BEFORE that instance's own declaration in the same
// scope. Declaration is interleaved with (and ordered by) source
// position, so at the point this assign is processed the target has
// no declared term yet. This -- like driving a child's internal net
// from outside its own scope at all -- isn't real synthesizable RTL
// (module ports are the only sanctioned cross-instance wiring
// mechanism), so it must throw rather than silently dropping the
// write and leaving the target fully unconstrained.
module sub;
  logic [1:0] x;
endmodule

module hier_continuous_assign_forward_ref;
  assign u1.x[0] = 1'b1;
  sub u1 ();
endmodule
