using System;
using System.IO;
using System.IO.Compression;
using System.Text;

namespace Ovation.Pipeline.FastqProcessor
{
    public class FastqWriter : ISequenceWriter
    {
        protected readonly StreamWriter streamWriter;

        private readonly FileStream outputStream;

        private readonly GZipStream gzipStream;

        protected readonly BufferedStream bufferedStream;

        protected readonly int bufferSize = 128 * 1024;

        private bool disposedValue;

        protected ulong sequencesWritten = 0;

        public ulong SequencesWritten => sequencesWritten;

        public FastqWriter(string fastq, bool gzipped)
        {
            var fileStreamOptions = new FileStreamOptions()
            {
                Mode = FileMode.Create,
                Access = FileAccess.ReadWrite,
                Share = FileShare.None,
                BufferSize = bufferSize,
            };

            if (gzipped == true)
            {
                outputStream = File.Open(fastq, fileStreamOptions);
                gzipStream = new GZipStream(outputStream, CompressionMode.Compress);
                bufferedStream = new BufferedStream(gzipStream, bufferSize);
            }
            else
            {
                outputStream = File.Open(fastq, fileStreamOptions);
                bufferedStream = new BufferedStream(outputStream, bufferSize);
            }

            streamWriter = new StreamWriter(bufferedStream, Encoding.ASCII, bufferSize, false);
        }

        public void WriteSequence(Sequence sequence)
        {
            try
            {
                streamWriter.WriteLine(sequence.Identifier);
                streamWriter.WriteLine(sequence.Read);
                streamWriter.WriteLine(sequence.Blank);
                streamWriter.WriteLine(sequence.Quality);

                sequencesWritten++;
            }
            catch (EndOfStreamException)
            {
                Console.Error.WriteLine("End of stream");
            }
        }

        protected void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    streamWriter?.Dispose();
                    bufferedStream?.Dispose();
                    gzipStream?.Dispose();
                    outputStream?.Dispose();
                }

                disposedValue = true;
            }
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
    }
}
