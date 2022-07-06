using System;

namespace Ovation.Pipeline.FastqProcessor
{
    public interface ISequenceReader : IDisposable
    {
        ulong SequencesRead { get; }

        bool ReadSequence(out Sequence sequence);

        double ApproximateCompletion { get; }
    }
}
