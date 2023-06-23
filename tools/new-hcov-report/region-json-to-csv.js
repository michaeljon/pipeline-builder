regions = require('./hcov-regions.json');

console.log('organism,sequence,gene,start,stop');

Object.keys(regions).forEach((rk, ri) => {
  region = regions[rk];
  sequence = region['sequence'];
  genes = region['genes'];

  Object.keys(genes).forEach((gk, gi) => {
    gene = genes[gk];

    console.log(`${rk},${sequence},${gk},${gene['start']},${gene['stop']}`);
  });
});
