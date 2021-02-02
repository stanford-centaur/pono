// example from Model Checking of Verilog RTL Using IC3 with Syntax-guided Abstraction
//              Aman Goel, Karem Sakallah

module uvw_example(clk, rst);
   parameter WIDTH=3;

   input clk, rst;
   reg [WIDTH-1:0]   u, v, w;

   wire [WIDTH-1:0]  one;
   assign one = 1;

   always @(posedge clk) begin
      if (rst) begin
         u <= one;
         v <= one;
         w <= one;
      end
      else begin
         u <= ((u < v) || (v < w)) ? u + v : v + one;
         v <= v + one;
         w <= w + one;
      end
   end

   always @* begin
      prop: assert((u + v) != one);
   end
endmodule // uv_example

