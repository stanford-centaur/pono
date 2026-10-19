// Concatenation-target LHS on a *procedural* (non-blocking) assignment
// (`{hi, lo} <= {a, b};`), the always_ff counterpart to the
// continuous-assign concat_lhs.sv test. begin_write() (shared by
// blocking/non-blocking assignment and ++/--) special-cases a
// top-level concatenation-target LHS the same way
// process_continuous_assign() does for the continuous-assign form.
module concat_lhs_next_state (input logic clk,
                              input logic rst,
                              input logic [3:0] a,
                              input logic [3:0] b);

  logic [3:0] hi, lo;

  always_ff @(posedge clk) begin
    if (rst) begin
      hi <= '0;
      lo <= '0;
    end else begin
      { hi, lo } <= { a, b };
    end
  end

  // {hi, lo} = {4'd10, 4'd11} = 8'd171 is reachable one cycle after
  // reset releases, once {hi, lo} tracks free inputs a/b.
  assert property (@(posedge clk) { hi, lo } != 8'd171);

endmodule
