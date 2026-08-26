// Integration test #2: a small single-write-port register file backed
// by a *packed* multi-dimensional array (`mem`) rather than a true SV
// unpacked array -- see Gap_UnpackedRegfileMemory in
// test_systemverilog_unsupported.cpp for why the latter isn't
// supported.  Combines: parameters, $clog2, a packed array-of-vectors
// read/written through a *runtime* address, and a one-cycle-delayed
// read-after-write invariant.
//
// Safety: if address X was written last cycle, reading address X this
// cycle must return the value that was written (single write port, so
// "last write wins" is well-defined).  Checked first with Bmc (can
// only report UNKNOWN, since it doesn't attempt an inductive proof),
// then with KInduction, which proves it (this invariant is
// 1-inductive).
module regfile_ctrl #(
    parameter int WIDTH = 8,
    parameter int DEPTH = 4
) (
    input  logic             clk,
    input  logic             rst,
    input  logic             we,
    input  logic [      1:0] waddr,
    input  logic [WIDTH-1:0] wdata,
    input  logic [      1:0] raddr,
    output logic [WIDTH-1:0] rdata
);

  localparam int IDX_W = $clog2(DEPTH);

  logic [DEPTH-1:0][WIDTH-1:0] mem;

  always_ff @(posedge clk) begin
    if (rst) mem <= '0;
    else if (we) mem[waddr] <= wdata;
  end

  assign rdata = mem[raddr];

  logic [1:0] waddr_prev;
  logic [WIDTH-1:0] wdata_prev;
  logic we_prev;

  always_ff @(posedge clk) begin
    if (rst) begin
      waddr_prev <= 2'd0;
      wdata_prev <= '0;
      we_prev    <= 1'b0;
    end else begin
      waddr_prev <= waddr;
      wdata_prev <= wdata;
      we_prev    <= we;
    end
  end

  assert property (@(posedge clk)
                    !(we_prev && raddr == waddr_prev) || rdata == wdata_prev);

endmodule
