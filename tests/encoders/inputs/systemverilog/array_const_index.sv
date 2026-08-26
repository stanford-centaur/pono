// Packed multi-dim array used as a 4-stage shift register, read/written
// through *constant* element indices: mem[0] free-runs like a plain
// counter, and mem[1..3] each shift in their left neighbor one cycle
// later.  Exercises ElementSelect on both the LHS (write) and RHS
// (read) of a packed array-of-vectors declaration, composed with
// always_ff sequencing (not a bare combinational one-liner).
//
// mem[0](k) = k-1 for k >= 1 (plain incrementer).  Each further stage
// lags one extra cycle behind the one before it (it reads the
// *pre-update* value of its neighbor), so mem[3](k) = k-4 for k >= 4:
// mem[3] is 0 for k <= 4 and first reaches 3 at cycle 7.
module array_const_index (input logic clk, input logic rst);

  logic [3:0][3:0] mem;

  always_ff @(posedge clk) begin
    if (rst) begin
      mem <= '0;
    end else begin
      mem[0] <= mem[0] + 4'd1;
      mem[1] <= mem[0];
      mem[2] <= mem[1];
      mem[3] <= mem[2];
    end
  end

  assert property (@(posedge clk) mem[3] != 4'd3);

endmodule
