// Submodule used by filelist_top.sv. Named by filelist.f rather than
// appearing directly in filelist_top.sv, to exercise --sv-filelist.
module child (
    input  logic [3:0] a,
    input  logic [3:0] b,
    output logic [3:0] sum
);
  assign sum = a + b;
endmodule
