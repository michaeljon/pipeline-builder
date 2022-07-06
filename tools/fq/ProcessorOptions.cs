using System;
using CommandLine;

namespace Ovation.Pipeline.FastqProcessor
{
    public abstract class ProcessorOptions
    {
        [Option('v', "verbose", Required = false, SetName = "Verbose", HelpText = "Set output to verbose messages.")]
        public bool Verbose { get; set; }

        [Option('d', "debug", Required = false, SetName = "Verbose", HelpText = "Show diagnostic output.  Can only use with --verbose.")]
        public bool Debug { get; set; }

        [Option("in1", Required = true, HelpText = "Forward strand (R1) input file")]
        public string ForwardInput { get; set; }

        [Option("in2", Required = true, HelpText = "Reverse strand (R2) input file")]
        public string ReverseInput { get; set; }

        [Option('f', "format", Required = true, HelpText = "Type of input file.")]
        public ReaderType Format { get; set; }

        [Option('l', "read-limit", Required = false, HelpText = "Limit the number of reads processed")]
        public ulong ReadLimit { get; set; } = ulong.MaxValue;
    }

    [Verb("split", HelpText = "Split the input R1/R2 files into partitions")]
    public class SplitOptions : ProcessorOptions
    {
        [Option("splits", Default = 128, HelpText = "Number of partitions to generate")]
        public int Splits { get; set; }

        [Option("output", Required = true, HelpText = "Output folder for R1/R2 partitions")]
        public string OutputFolder { get; set; }
    }

    [Verb("sample", HelpText = "Select a sample of the input reads into a new R1/R2 pair")]
    public class SampleOptions : ProcessorOptions
    {
        [Option("sample-rate", Default = 128, HelpText = "Select one from every sample-rate sequences to generate output files")]
        public int SampleRate { get; set; }

        [Option("out1", Required = false, HelpText = "Forward strand (R1) output file. Will default to filename_sampled.")]
        public string ForwardOutput { get; set; }

        [Option("out2", Required = false, HelpText = "Reverse strand (R2) output file. Will default to filename_sampled.")]
        public string ReverseOutput { get; set; }
    }
}
