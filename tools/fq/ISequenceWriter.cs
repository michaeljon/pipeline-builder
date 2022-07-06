using System;

namespace Ovation.Pipeline.FastqProcessor
{
    public interface ISequenceWriter : IDisposable
    {
        ulong SequencesWritten { get; }

        void WriteSequence(Sequence sequence);
    }
}
