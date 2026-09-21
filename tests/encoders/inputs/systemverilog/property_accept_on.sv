// `accept_on` as a top-level property connective. LRM 1800-2009
// F.5.3.1: it holds if the operand does, or if some cycle has the
// abort condition and the operand holds on the word truncated just
// before it and padded with true. For an operand that reduces to a
// per-cycle check, that padding cannot rescue a failure already
// observed and forgives everything from the abort onward -- so the
// check is the operand's own, waived across the cycles an attempt
// spans if the condition held in any of them.
//
// Here `abrt` is free and can stay low, so the waiver need not apply
// and the free `a` refutes it immediately.
module property_accept_on (
    input logic clk,
    input logic abrt,
    input logic a
);

  assert property (@(posedge clk) accept_on (abrt) a);

endmodule
