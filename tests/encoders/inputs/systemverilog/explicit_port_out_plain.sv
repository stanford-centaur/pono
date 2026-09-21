// explicit_port_out.sv with an ordinary output port, which is the
// control: whatever the explicit form does, it has to do this.
module explicit_port_out_plain_sub (
    input logic clk,
    input logic d,
    output logic o
);
  always_ff @(posedge clk) o <= d;
endmodule

module explicit_port_out_plain (
    input logic clk,
    input logic d
);

  logic q;
  logic started;
  initial started = 1'b0;
  always_ff @(posedge clk) started <= 1'b1;

  explicit_port_out_plain_sub u (
      .clk(clk),
      .d  (d),
      .o  (q)
  );

  assert property (@(posedge clk) started |-> q == $past(d));

endmodule
