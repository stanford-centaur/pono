// A recursive function whose depth is a runtime value. Not
// synthesizable -- it describes no finite piece of hardware, and
// synthesis tools reject it too -- so it is out of scope rather than
// unimplemented. (Recursion to a *constant* depth is a different
// question: that elaborates statically and could be unrolled.)
module unbounded_recursion (input logic clk, input logic [7:0] d,
                            output logic [7:0] q);

  function automatic [7:0] f(input [7:0] n);
    if (n == 8'd0) f = 8'd0;
    else f = n + f(n - 8'd1);
  endfunction

  assign q = f(d);

  assert property (@(posedge clk) q >= 8'd0);

endmodule
