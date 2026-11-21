// An `interface` with `modport`s, accessed through a modport-qualified
// port (`bus_if.master b`) rather than the plain signal-bundle style
// in interface_bundle.sv. The `task` is declared but unused by the
// design under test; it is only present to show a task inside an
// interface doesn't break elaboration.
interface bus_if;
  logic [3:0] data;
  modport master(output data);
  modport slave(input data);
  task automatic reset_bus();
    data = 4'd0;
  endtask
endinterface

module mst (bus_if.master b, input logic clk, input logic rst);
  always_ff @(posedge clk) begin
    if (rst) b.data <= 4'd0;
    else b.data <= b.data + 4'd1;
  end
endmodule

module interface_modport_task (input logic clk, input logic rst);
  bus_if b ();
  mst m (
      b.master,
      clk,
      rst
  );
  assert property (@(posedge clk) b.data != 4'd3);
endmodule
