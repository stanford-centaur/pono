// Casting a plain integer expression to an enum type (`state_t'(x)`)
// -- an explicit Conversion into an integral-but-enumerated type,
// combined with a register driven by the cast result.  `sel` is a
// free 2-bit input; casting it directly into `state_t` and registering
// it means BMC can reach the out-of-declared-range encoding 2'b11 (not
// one of IDLE/REQ/ACK) immediately at cycle 1, since 2'b11 has no
// symbolic-integral restriction once cast.
typedef enum logic [1:0] {
  IDLE = 2'b00,
  REQ  = 2'b01,
  ACK  = 2'b10
} state_t;

module enum_cast (input logic clk, input logic rst, input logic [1:0] sel);

  state_t st;

  always_ff @(posedge clk) begin
    if (rst) st <= IDLE;
    else st <= state_t'(sel);
  end

  assert property (@(posedge clk) st != 2'b11);

endmodule
