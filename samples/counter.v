module counter(clk, in, out);
   parameter WIDTH=8;

   input                    clk, in;
   output reg [WIDTH-1:0]   out;

   always @(posedge clk) begin
      out <= out + in;
   end

   reg initstate = 1'b1;

   always @(posedge clk) begin
      initstate <= 1'b0;
   end

   always @* begin
      // unsafe
      assert(out < 10);
   end

endmodule // counter
