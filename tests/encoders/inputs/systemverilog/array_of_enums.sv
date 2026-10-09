// A packed array whose element type is a user-defined enum (rather
// than a plain bit-vector), read/written through a constant index.
// Exercises the combination of packed-array ElementSelect with an
// enum base type.
//
// slot[0] cycles IDLE -> REQ -> ACK -> IDLE -> ...; it first reaches
// ACK at cycle 3 (one reset cycle + two transitions).
typedef enum logic [1:0] {
  IDLE = 2'b00,
  REQ  = 2'b01,
  ACK  = 2'b10
} state_t;

module array_of_enums (input logic clk, input logic rst);

  state_t [1:0] slot;

  always_ff @(posedge clk) begin
    if (rst) begin
      slot[0] <= IDLE;
      slot[1] <= IDLE;
    end else begin
      case (slot[0])
        IDLE: slot[0] <= REQ;
        REQ: slot[0] <= ACK;
        ACK: slot[0] <= IDLE;
        default: slot[0] <= IDLE;
      endcase
      slot[1] <= slot[0];
    end
  end

  assert property (@(posedge clk) slot[0] != ACK);

endmodule
