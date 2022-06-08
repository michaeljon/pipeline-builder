const fs = require('fs');
const f = fs.readFileSync('./human/hg38-chromosomeSizes.json');
const d = JSON.parse(f);

const neucleotides = ['A', 'C', 'T', 'G'];
const ncount = neucleotides.length;

// https://academic.oup.com/hmg/article/19/R2/R131/640505
const variantChance = 0.999;
const indelChance = 0.18;

const randomRead = () => neucleotides[parseInt(Math.random() * ncount)];

function generateIndel() {
  const l = parseInt(Math.random() * (16 - 2) + 2, 10);
  const indel = Array(l);
  for (let n = 0; n < l; n++) {
    indel[n] = randomRead();
  }
  return indel.join('');
}

const writeVariant = (chr, locus, r, a) => console.log(`chr${chr}\t${locus}\t.\t${r}\t${a}\t0.0\tDP=0`);

let variantCount = 0;
let indelCount = 0;

function generateVariants(chr, size) {
  for (let locus = 1; locus <= size; locus++) {
    if (Math.random() > variantChance) {
      variantCount++;

      if (Math.random() > indelChance) {
        // is snp
        writeVariant(chr, locus, randomRead(), randomRead());
      } else {
        // is indel
        const indel = generateIndel();

        indelCount++;
        if (Math.random() > 0.5) {
          // insertion
          writeVariant(chr, locus, indel.charAt(0), indel);
        } else {
          // deletion
          writeVariant(chr, locus, indel, indel.charAt(0));
        }

        locus += indel.length;
      }
    }
  }
}

Object.keys(d).forEach(k => generateVariants(k, parseInt(d[k], 10)));

console.log(`Number of variants ${variantCount}`);
console.log(`Number of indels ${indelCount}`);
