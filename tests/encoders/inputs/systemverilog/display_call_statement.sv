// Every class of system task slang registers that only reports:
// console and file output, severity, simulation control, waveform
// dumping, time formatting. None of them touches design state, so
// the register they surround is still proved to follow its input.
//
// The tasks that do touch state -- $readmemh, $deposit, $cast,
// $sformat and the rest that take an argument by reference -- are
// refused instead; see system_task_output_arg.sv.
module display_call_statement (
    input logic       clk,
    input logic       rst,
    input logic [3:0] a
);

  logic [3:0] b;
  logic [7:0] trace[0:3];

  always_ff @(posedge clk) begin
    if (rst) b <= 4'd0;
    else begin
      // Console output, in every radix and both flavours.
      $display("a = %0d", a);
      $displayb(a);
      $displayo(a);
      $displayh(a);
      $write("a = %0d\n", a);
      $writeb(a);
      $writeo(a);
      $writeh(a);

      // The end-of-timestep and value-change variants.
      $strobe("a = %0d", a);
      $strobeb(a);
      $monitor("a = %0d", a);
      $monitorh(a);
      $monitoron;
      $monitoroff;

      // The same families writing to a file descriptor.
      $fdisplay(32'h8000_0001, "a = %0d", a);
      $fwrite(32'h8000_0001, "a = %0d", a);
      $fstrobe(32'h8000_0001, "a = %0d", a);
      $fmonitor(32'h8000_0001, "a = %0d", a);

      // Severity, which reports without changing anything.
      $info("just so you know");
      $warning("a is %0d", a);
      $error("a is %0d", a);
      $fatal(0, "a is %0d", a);

      // Simulation control: there is no run to start or stop here.
      $finish;
      $stop;
      $exit;

      // Waveform dumping, which only observes.
      $dumpfile("dump.vcd");
      $dumpvars(0, display_call_statement);
      $dumpon;
      $dumpoff;
      $dumpall;
      $dumplimit(1024);
      $dumpflush;

      // Time formatting and the shell escape.
      $timeformat(-9, 2, " ns", 10);
      $printtimescale;
      $system("true");

      // Writing a memory *out* only reads it.
      $writememh("trace.hex", trace);
      $writememb("trace.bin", trace);

      b <= a;
    end
  end

  assert property (@(posedge clk) rst || $past(rst) || b == $past(a));

endmodule
