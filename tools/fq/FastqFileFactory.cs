using System;

namespace Ovation.Pipeline.FastqProcessor
{
    public static class FastqFileFactory
    {
        public static ISequenceReader CreateReader(ReaderType format, string fastq)
        {
            return format switch
            {
                ReaderType.Fastq => new FastqReader(fastq, false),
                ReaderType.FastqGz => new FastqReader(fastq, true),

                _ => throw new InvalidOperationException($"could not determine file type of {fastq}")
            };
        }

        public static ISequenceWriter CreateWriter(ReaderType format, string fastq)
        {
            return format switch
            {
                ReaderType.Fastq => new FastqWriter(fastq, false),
                ReaderType.FastqGz => new FastqWriter(fastq, true),

                _ => throw new InvalidOperationException($"could not determine file type of {fastq}")
            };
        }
    }
}
