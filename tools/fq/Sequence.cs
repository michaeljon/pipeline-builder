using System;
using System.Text;

namespace Ovation.Pipeline.FastqProcessor
{
    public class Sequence
    {
        public string Identifier { get; }

        public string Read { get; }

        public string Blank { get; }

        public string Quality { get; }

        public Sequence(string identifer, string read, string blank, string quality)
        {
            Identifier = identifer;
            Read = read;
            Blank = blank;
            Quality = quality;
        }

        public override string ToString()
        {
            var sb = new StringBuilder("sequence: \n");

            sb.AppendLine(Identifier);
            sb.AppendLine(Read);
            sb.AppendLine(Blank);
            sb.AppendLine(Quality);

            return sb.ToString();
        }
    }
}
