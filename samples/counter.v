module counter(clk, in, out);
   parameter WIDTH=8;

   input                    clk, in;
   output wire [WIDTH-1:0]   out;

   reg [WIDTH:0]             internal = 0;
   assign out = internal;

   always @(posedge clk) begin
      internal <= internal + in;
   end

   reg initstate = 1'b1;

   always @(posedge clk) begin
      initstate <= 1'b0;
   end

   always @* begin
      // unsafe
      assert(internal < 10);
   end

endmodule // counter
