// A user-defined packed struct type used as a *module port* type
// (not just an internal variable), read through two of its fields in
// the same cycle.  Both fields are free (4-bit) inputs, so their sum
// can reach at most 15+15=30.  Cycle 1 is the reset value (acc == 0,
// din not yet sampled); the earliest transition that can actually use
// din is the one into cycle 2, so that is the earliest falsifying
// bound.
typedef struct packed {
  logic [3:0] a;
  logic [3:0] b;
} sum_in_t;

module typedef_struct_port (
    input logic clk,
    input logic rst,
    input sum_in_t din
);

  logic [4:0] acc;

  always_ff @(posedge clk) begin
    if (rst) acc <= 5'd0;
    else acc <= {1'b0, din.a} + {1'b0, din.b};
  end

  assert property (@(posedge clk) acc < 5'd30);

endmodule
