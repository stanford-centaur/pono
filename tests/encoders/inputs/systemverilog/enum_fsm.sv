// Packed enum used as an always_ff state register, driving a small
// FSM (IDLE -> REQ -> ACK -> IDLE) through a `case` on the enum
// value.  `go` is a free input; BMC can release it to reach ACK.
// Paired with enum_fsm_holds.sv, which ties `go` low forever so the
// FSM can never leave IDLE and the same property genuinely holds.
// Deliberately has no `default:` case arm (2'b11 is an unreachable
// encoding here, so the register simply holds its value in that
// state) -- a separate, distinct encoder bug makes a `case`
// statement's `default:` arm apply unconditionally alongside
// whichever other item already matched, rather than only when none
// of them did, which would otherwise stick this FSM at IDLE forever.
typedef enum logic [1:0] {
  IDLE = 2'b00,
  REQ  = 2'b01,
  ACK  = 2'b10
} state_t;

module enum_fsm (input logic clk, input logic rst, input logic go);

  state_t st;

  always_ff @(posedge clk) begin
    if (rst) begin
      st <= IDLE;
    end else begin
      case (st)
        IDLE: st <= go ? REQ : IDLE;
        REQ: st <= ACK;
        ACK: st <= IDLE;
      endcase
    end
  end

  assert property (@(posedge clk) st != ACK);

endmodule
