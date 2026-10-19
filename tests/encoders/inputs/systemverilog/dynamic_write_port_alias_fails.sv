// The paired refutation for dynamic_write_port_alias.sv: when the
// index selects the low half, the high half is not the one written,
// so claiming it was must be refuted. Without this, a write that
// reached neither half would satisfy the holds case.
module dynamic_write_port_alias_fails (
    input logic clk,
    input logic rst,
    input logic [2:0] i
);

  logic [3:0] hi, lo;

  alias_driver2 u (
      .clk(clk),
      .rst(rst),
      .i  (i),
      .p  ({hi, lo})
  );

  assert property (@(posedge clk) (!rst && !$past(rst) && $past(i) < 3'd4)
                   |-> hi[$past(i)] == 1'b1);

endmodule

module alias_driver2 (
    input logic clk,
    input logic rst,
    input logic [2:0] i,
    output logic [7:0] p
);
  always_ff @(posedge clk) begin
    if (rst) p <= 8'd0;
    else p[i] <= 1'b1;
  end
endmodule
