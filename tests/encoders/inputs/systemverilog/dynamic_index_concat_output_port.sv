// A submodule's output port connected to a concatenation
// (`.out({a, bus[idx]})`) where one operand resolves fine but the
// other is a dynamic (runtime-variable) bit-select. The
// concatenation-target output-port path must throw as soon as any
// operand's resolve_lvalue() fails, rather than silently dropping the
// whole multi-piece write and leaving both `a` and `bus` unconstrained.
module child (output logic [1:0] out, input logic [1:0] in);

  assign out = in;

endmodule

module dynamic_index_concat_output_port (input logic [1:0] idx,
                                         input logic [1:0] in);

  logic a;
  logic [3:0] bus;

  child c (.out({ a, bus[idx] }), .in(in));

endmodule
