module mem(clk, wr, addr, data_in, data_out);
   parameter WIDTH      = 8;
   parameter DEPTH      = 16;
   parameter ADDR_WIDTH = $clog2(DEPTH);

   input                   clk, wr;
   input  [ADDR_WIDTH-1:0] addr;
   input  [WIDTH-1:0]      data_in;
   output [WIDTH-1:0]      data_out;

   reg [WIDTH-1:0]         mem [DEPTH-1:0];

   always @(posedge clk) begin
      if (wr) begin
         mem[addr] <= data_in;
      end
   end

   assign data_out = mem[addr];
endmodule // mem
