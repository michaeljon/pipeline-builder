const chromosomeSizes = require('./GRCh38-p14-chromosomeSizes.json');

const computeIntervals = () => {
  const segmentSize = 50000000;
  const overlapSize = 50000;
  const lastBlockMax = 50000000 * 0.25;

  const intervals = [];
  chromosomeSizes.forEach(c => {
    let remainder = c.length - segmentSize;
    const segments = Math.ceil(c.length / segmentSize);
    let segment = 0;

    while (remainder > lastBlockMax) {
      const lower = segment * segmentSize + 1;
      const upper = (segment + 1) * segmentSize;

      if (lower > 1) {
        intervals.push(`${c.accession}:${lower - overlapSize - 1}-${lower + overlapSize}`);
      }

      intervals.push(`${c.accession}:${lower}-${upper}`);

      segment++;
      remainder -= segmentSize;
    }

    let lower = (segments - 2) * segmentSize;

    if (remainder <= 0) {
      lower = (segments - 1) * segmentSize;
    }

    if (lower % 10 == 0) lower += 1;

    if (lower > 1) {
      intervals.push(`${c.accession}:${lower - overlapSize - 1}-${lower + overlapSize}`);
    }

    intervals.push(`${c.accession}:${lower}-${c.length}`);
  });

  return intervals;
};

console.log(computeIntervals());
