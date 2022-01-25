const chromosomeSizes = require('./chromosomeSizes.json');

const computeIntervals = () => {
  const segmentSize = 50000000;
  const lastBlockMax = 50000000 * 0.25;

  const intervals = [];
  Object.keys(chromosomeSizes).forEach((c) => {
    let remainder = chromosomeSizes[c] - segmentSize;
    const segments = Math.ceil(chromosomeSizes[c] / segmentSize);
    let segment = 0;

    while (remainder > lastBlockMax) {
      const lower = segment * segmentSize + 1;
      const upper = (segment + 1) * segmentSize;

      intervals.push(`"chr${c}:${lower}-${upper}"`);

      segment++;
      remainder -= segmentSize;
    }

    if (remainder > 0) {
      intervals.push(`"chr${c}:${(segments - 2) * segmentSize}-${chromosomeSizes[c]}"`);
    } else {
      intervals.push(`"chr${c}:${(segments - 1) * segmentSize}-${chromosomeSizes[c]}"`);
    }
  });

  return intervals;
};

console.log(computeIntervals());
