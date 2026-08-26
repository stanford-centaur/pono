// `a[+]` (`[*1:$]`, unbounded consecutive repetition), same
// architectural boundary as unbounded_repeat_star.sv.
module unbounded_repeat_plus (input logic clk, input logic rst, input logic a);

  assert property (@(posedge clk) a[+] |-> 1'b0);

endmodule
