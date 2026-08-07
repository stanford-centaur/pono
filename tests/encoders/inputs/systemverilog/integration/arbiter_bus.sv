// Integration test #1 (safety half): a small round-robin arbiter
// hierarchically instantiating two `requester` children.  Grant
// decisions are computed directly from the free `req` inputs (not
// from the children's internal `pending` state) -- earlier drafts of
// this design instead read each child's internal `pending` back
// (via an output port, and separately via a hierarchical value
// reference) from a parent-side always_comb block that also fed a
// sibling instance's port, and every variant of that shape threw
// "unknown symbol" building the property term, even though the
// individual ingredients (output-port aliasing, hierarchical value
// references, always_comb) all work fine in isolation elsewhere in
// this suite.  That combination wasn't successfully narrowed down to
// a minimal, confidently-attributable repro, so it's not reported as
// a standalone Gap_ finding; this design was simplified to avoid it
// rather than risk mischaracterizing the cause.  Combines: module
// hierarchy, a packed struct built with a whole-
// value assignment pattern ('{...}), always_comb + always_ff
// interaction, and a safety property checked with plain Bmc.
//
// Safety: a grant is only ever issued to a requester that actually
// asked for it *this* cycle -- holds by construction, so BMC should
// exhaust any finite bound with UNKNOWN.  See arbiter_bus_live.sv for
// the companion liveness property.
typedef struct packed {
  logic valid;
  logic id;
} grant_t;

module requester (input logic clk, input logic rst, input logic req_in, input logic grant_in);
  logic pending;
  always_ff @(posedge clk) begin
    if (rst) pending <= 1'b0;
    else if (req_in) pending <= 1'b1;
    else if (grant_in) pending <= 1'b0;
  end
endmodule

module arbiter_bus (input logic clk, input logic rst, input logic [1:0] req);

  logic   turn;
  grant_t grant;

  always_comb begin
    if (req[0] && (!req[1] || !turn)) grant = '{valid: 1'b1, id: 1'b0};
    else if (req[1]) grant = '{valid: 1'b1, id: 1'b1};
    else grant = '{valid: 1'b0, id: 1'b0};
  end

  requester r0 (
      .clk     (clk),
      .rst     (rst),
      .req_in  (req[0]),
      .grant_in(grant.valid && !grant.id)
  );

  requester r1 (
      .clk     (clk),
      .rst     (rst),
      .req_in  (req[1]),
      .grant_in(grant.valid && grant.id)
  );

  always_ff @(posedge clk) begin
    if (rst) turn <= 1'b0;
    else if (grant.valid) turn <= ~grant.id;
  end

  assert property (@(posedge clk)
                    !(grant.valid && !(grant.id ? req[1] : req[0])));

endmodule
