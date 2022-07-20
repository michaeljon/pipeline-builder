using System;
using System.IO;
using System.IO.Compression;
using System.Text;

namespace Ovation.Pipeline.FastqProcessor
{
    public class FastqReader : ISequenceReader
    {
        protected readonly StreamReader streamReader;

        private readonly FileStream inputStream;

        private readonly GZipStream gzipStream;

        protected readonly BufferedStream bufferedStream;

        protected readonly int bufferSize = 128 * 1024;

        private bool disposedValue;

        protected ulong sequencesRead = 0;

        public ulong SequencesRead => sequencesRead;

        public double ApproximateCompletion =>
            100.0 * inputStream.Position / inputStream.Length;

        public FastqReader(string fastq, bool gzipped)
        {
            var fileStreamOptions = new FileStreamOptions()
            {
                Mode = FileMode.Open,
                BufferSize = bufferSize,
            };

            if (gzipped == true)
            {
                inputStream = File.Open(fastq, fileStreamOptions);
                gzipStream = new GZipStream(inputStream, CompressionMode.Decompress);
                bufferedStream = new BufferedStream(gzipStream, bufferSize);
            }
            else
            {
                inputStream = File.Open(fastq, fileStreamOptions);
                bufferedStream = new BufferedStream(inputStream, bufferSize);
            }

            streamReader = new StreamReader(bufferedStream, Encoding.ASCII, false, bufferSize);
        }

        public bool ReadSequence(out Sequence sequence)
        {
            try
            {
                if (streamReader.EndOfStream == true)
                {
                    Console.Error.WriteLine("End of stream");
                    sequence = null;
                    return false;
                }

                var identifier = streamReader.ReadLine();
                var read = streamReader.ReadLine();
                var blank = streamReader.ReadLine();
                var quality = streamReader.ReadLine();

                sequence = new Sequence(identifier, read, blank, quality);
                sequencesRead++;
                return true;
            }
            catch (EndOfStreamException)
            {
                Console.Error.WriteLine("End of stream");
                sequence = null;
                return false;
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine(ex.Message);
                sequence = null;
                Environment.Exit(1);

                // unreachable code
                return false;
            }
        }

        protected void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    streamReader?.Dispose();
                    bufferedStream?.Dispose();
                    gzipStream?.Dispose();
                    inputStream?.Dispose();
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
