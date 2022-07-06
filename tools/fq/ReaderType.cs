using System;

namespace Ovation.Pipeline.FastqProcessor
{
    [Flags]
    public enum ReaderType
    {
        Fastq = 1,

        FastqGz = 2,

        UnalignedReaders = Fastq | FastqGz,

        AllReaders = UnalignedReaders
    }
}