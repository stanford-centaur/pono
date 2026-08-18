// Same FSM shape as enum_fsm.sv (including the missing `default:`
// arm, for the same reason -- see that file's comment), but `go` is
// hardwired low so the state machine can never leave IDLE.  The
// property `st != ACK` genuinely holds forever here -- BMC must
// exhaust the bound with UNKNOWN rather than manufacture a
// counterexample, which is exactly the failure mode a test that only
// ever falsifies its property could never catch.
typedef enum logic [1:0] {
  IDLE = 2'b00,
  REQ  = 2'b01,
  ACK  = 2'b10
} state_t;

module enum_fsm_holds (input logic clk, input logic rst);

  state_t st;
  logic go;

  assign go = 1'b0;

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
