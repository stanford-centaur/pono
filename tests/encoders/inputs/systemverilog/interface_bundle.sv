// A plain signal-bundle `interface` (no modports, no tasks) wiring a
// producer module to the top level -- the most common "just group
// some wires" use of `interface` in synthesizable RTL, distinct from
// the modport/task-bearing interfaces covered as out-of-scope in
// test_systemverilog_unsupported.cpp.  bus.data(k) = k-1, so it first
// equals 3 at cycle 4, *if* interface instances/ports are walked like
// ordinary module instances.
interface simple_bus;
  logic [3:0] data;
  logic       valid;
endinterface

module bus_producer (simple_bus bus, input logic clk, input logic rst);
  always_ff @(posedge clk) begin
    if (rst) begin
      bus.data  <= 4'd0;
      bus.valid <= 1'b0;
    end else begin
      bus.data  <= bus.data + 4'd1;
      bus.valid <= 1'b1;
    end
  end
endmodule

module interface_bundle (input logic clk, input logic rst);
  simple_bus bus ();
  bus_producer p (
      bus,
      clk,
      rst
  );
  assert property (@(posedge clk) bus.data != 4'd3);
endmodule
