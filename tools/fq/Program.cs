using System;
using System.IO;
using System.Threading.Tasks;
using CommandLine;

namespace Ovation.Pipeline.FastqProcessor
{
    internal class Program
    {
        static void Main(string[] args)
        {
            Parser parser = new(config =>
                {
                    config.AutoHelp = true;
                    config.AutoVersion = true;
                    config.CaseInsensitiveEnumValues = true;
                    config.EnableDashDash = false;
                    config.HelpWriter = Console.Error;
                }
            );

            parser.ParseArguments<SplitOptions, SampleOptions>(args)
                .WithParsed<SplitOptions>((SplitOptions o) => RunSplitsDriver(o))
                .WithParsed<SampleOptions>((SampleOptions o) => RunSamplesDriver(o))
                .WithNotParsed(errors => {});
        }

        private static int RunSplitsDriver(SplitOptions o)
        {
            var t1 = Task.Run(() => RunSplits(o.Format, o.ReadLimit, o.Splits, o.ForwardInput, o.OutputFolder));
            var t2 = Task.Run(() => RunSplits(o.Format, o.ReadLimit, o.Splits, o.ReverseInput, o.OutputFolder));

            Task.WaitAll(t1, t2);

            Console.WriteLine("{0} sequences processed", t1.Result + t2.Result);

            return 0;
        }

        private static ulong RunSplits(ReaderType readerType, ulong readLimit, int splits, string inputFile, string outputFolder)
        {
            using ISequenceReader sequenceReader = FastqFileFactory.CreateReader(readerType, inputFile);

            var basename = Path.GetFileName(inputFile)
                .Replace(".fastq.gz", "")
                .Replace(".fastq", "");

            var filenames = new string[splits];
            var writers = new ISequenceWriter[splits];
            for (var w = 0; w < splits; w++)
            {
                filenames[w] = Path.Combine(outputFolder, string.Format("{0}_{1:d3}.{2}", basename, w + 1, readerType == ReaderType.Fastq ? "fastq" : "fastq.gz"));
                writers[w] = FastqFileFactory.CreateWriter(readerType, filenames[w]);
            }

            var split = 0;
            while (sequenceReader.SequencesRead < readLimit && sequenceReader.ReadSequence(out Sequence sequence))
            {
                writers[split++ % splits].WriteSequence(sequence);
            }

            for (var w = 0; w < splits; w++)
            {
                Console.WriteLine("{0} sequences written to {1}", writers[w].SequencesWritten, filenames[w]);
                writers[w].Dispose();
            }

            return sequenceReader.SequencesRead;
        }

        private static int RunSamplesDriver(SampleOptions o)
        {
            var t1 = Task.Run(() => RunSamples(o.Format, o.ReadLimit, o.SampleRate, o.ForwardInput, o.ForwardOutput));
            var t2 = Task.Run(() => RunSamples(o.Format, o.ReadLimit, o.SampleRate, o.ReverseInput, o.ReverseOutput));

            Task.WaitAll(t1, t2);

            Console.WriteLine("{0} sequences processed", t1.Result + t2.Result);

            return 0;
        }

        private static ulong RunSamples(ReaderType readerType, ulong readLimit, int sampleRate, string inputFile, string outputFile)
        {
            var targetFile = outputFile;

            if (string.IsNullOrEmpty(targetFile))
            {
                var basename = Path.GetFileName(outputFile)
                    .Replace(".fastq.gz", "")
                    .Replace(".fastq", "");
                targetFile = string.Format("{0}_sampled.{1}", basename, readerType == ReaderType.Fastq ? "fastq" : "fastq.gz");
            }

            using ISequenceReader sequenceReader = FastqFileFactory.CreateReader(readerType, inputFile);
            using ISequenceWriter sequenceWriter = FastqFileFactory.CreateWriter(readerType, targetFile);

            var read = 0UL;
            var written = 0UL;
            while (sequenceReader.SequencesRead < readLimit && sequenceReader.ReadSequence(out Sequence sequence))
            {
                if (read++ % (ulong)sampleRate == 0)
                {
                    sequenceWriter.WriteSequence(sequence);
                    written++;
                }
            }

            Console.WriteLine("{0} sequences written to {1}", written, outputFile);

            return read;
        }
    }
}
