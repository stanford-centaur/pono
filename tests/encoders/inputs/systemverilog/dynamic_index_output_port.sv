// A submodule's output port connected to a dynamic (runtime-variable)
// bit-select of an internal signal (`.out(bus[idx])`). Unlike a
// procedural assignment, a port connection is a structural,
// elaboration-time binding with no per-cycle write / mux fallback, so
// resolve_lvalue()'s dynamic-index failure here must throw a clear
// PonoException rather than silently dropping the child's output
// write and leaving `bus` fully unconstrained.
module child (output logic out, input logic in);

  assign out = in;

endmodule

module dynamic_index_output_port (input logic [1:0] idx,
                                  input logic in);

  logic [3:0] bus;

  child c (.out(bus[idx]), .in(in));

endmodule
