// A recursive function whose depth is a runtime value. Not
// synthesizable -- there is no finite piece of hardware it
// describes, and synthesis tools reject it too -- so it is out of
// scope here rather than unimplemented. (A constant-depth recursion
// is a different question: that one elaborates statically and could
// be unrolled. This is the unbounded kind.)
//
// It currently reports the recursion as something else entirely (a
// local variable written under a runtime condition), because the
// walk reaches local_var_both_arms.sv's bug on the way down. That is
// tracked there; all this asserts is the refusal.
module unbounded_recursion (input logic clk, input logic [7:0] d,
                            output logic [7:0] q);

  function automatic [7:0] f(input [7:0] n);
    if (n == 8'd0) f = 8'd0;
    else f = n + f(n - 8'd1);
  endfunction

  assign q = f(d);

  assert property (@(posedge clk) q >= 8'd0);

endmodule
