// A recursive function. Recursion has no point at which inlining
// stops, so refusing it is right -- but it is refused as "the local
// variable 'f' is written only under a runtime condition", which
// names the function's own return variable and says nothing about
// the recursion. The message is the gap here, not the refusal.
//
// Intended: a thrown error that says what is actually wrong. The
// test matches on the word, so it goes green when the message names
// recursion and stays red for any wording that does not.
module recursive_subroutine (input logic clk, input logic [7:0] d,
                             output logic [7:0] q);

  function automatic [7:0] f(input [7:0] n);
    if (n == 8'd0) f = 8'd0;
    else f = n + f(n - 8'd1);
  endfunction

  assign q = f(d);

  assert property (@(posedge clk) q >= 8'd0);

endmodule
