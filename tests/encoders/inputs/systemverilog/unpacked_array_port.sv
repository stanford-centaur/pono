// An unpacked array crossing a module boundary, in both directions
// and at the top level.
//
// A whole array has no bits to splice across the port-connection
// machinery's alias segments, so each case is one term shared by
// both sides: an input port takes the parent's array term, and an
// output port makes the child's port variable and the parent's array
// the same state variable, which the child's writes then drive.
//
// `seen` reads a memory passed down into a submodule; `driven` is a
// memory written by one. unpacked_array_port_fails.sv is the paired
// refutation for the second, since an array nothing constrains would
// satisfy any expected value at all.
module unpacked_array_port (
    input logic clk,
    input logic rst,
    input logic [1:0] i,
    input logic [7:0] d,
    input logic [7:0] ext[0:3]
);

  logic [7:0] mem[0:3];
  logic [7:0] driven[0:3];
  logic [7:0] seen;
  logic [7:0] from_ext;

  always_ff @(posedge clk) begin
    if (rst) mem <= '{default: 0};
    else mem[i] <= d;
  end

  array_reader u_read (
      .clk(clk),
      .m  (mem),
      .i  (i),
      .o  (seen)
  );

  array_reader u_ext (
      .clk(clk),
      .m  (ext),
      .i  (i),
      .o  (from_ext)
  );

  array_writer u_write (
      .clk(clk),
      .rst(rst),
      .i  (i),
      .d  (d),
      .m  (driven)
  );

  assert property (@(posedge clk)
                   (seen == mem[i] && from_ext == ext[i]
                    && (!(!rst && !$past(rst) && $past(i) == i)
                        || driven[i] == $past(d))));

endmodule

module array_reader (
    input logic clk,
    input logic [7:0] m[0:3],
    input logic [1:0] i,
    output logic [7:0] o
);
  assign o = m[i];
endmodule

module array_writer (
    input logic clk,
    input logic rst,
    input logic [1:0] i,
    input logic [7:0] d,
    output logic [7:0] m[0:3]
);
  always_ff @(posedge clk) begin
    if (rst) m <= '{default: 0};
    else m[i] <= d;
  end
endmodule
