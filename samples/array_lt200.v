module array_pred(clk, wr_addr, rd_addr, data, zeros);

   input clk;
   input [7:0] wr_addr, rd_addr;
   input [31:0] data;
   output reg [31:0]  zeros;

   reg [31:0]   mem [255:0];

   always @(posedge clk) begin
      mem[wr_addr] <= data < 32'd200 ? data : mem[wr_addr];
      zeros <= 1'b0;
   end

   always @* begin
      lt_200: assert(mem[rd_addr] < 32'd200);
   end

endmodule // array_pred
