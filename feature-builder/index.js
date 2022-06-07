const fs = require('fs');


// list of primary regions in fasta
// grep -e '^>[^ ]* Homo sapiens chromosome [1-9XY]*, GRCh38.p14 Primary Assembly' GCF_000001405.40_GRCh38.p14_genomic.fna

// sanity-check regions in fasta
// / awk '$1~/>/{a[$1]++}END{for (e in a) { if (a[e] > 1) { print e, a[e] }}}' GCF_000001405.40_GRCh38.p14_genomic.fna

// count of features per region
// awk '!($1~/#/){a[$1]++}END{for (e in a) { print e, a[e] }}' GCF_000001405.40_GRCh38.p14_genomic.gff  | sort -k2 -g -r

// index int8range layout

//                             bit position

//   6 6 6 6 5 5 5 5 5 5 5 5 5 5 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3
//   3 2 1 0 9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2
// [ p|  reserved   |  chromosome   |           sequence             ]

//   3 3 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
//   1 0 9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2 1 0
// [                            locus                                ]


const FASTA =
  '/Users/michaeljon/Downloads/grch38p14 genome_assemblies_genome_gff/ncbi-genomes-2022-06-06/GCF_000001405.40_GRCh38.p14_genomic.fna';
const GFF =
  '/Users/michaeljon/Downloads/grch38p14 genome_assemblies_genome_gff/ncbi-genomes-2022-06-06/GCF_000001405.40_GRCh38.p14_genomic.gff';

let sequenceNumber = 1024;
let lineNo = 0;
const sequences = {};

function processFeature(line) {
  if (line.startsWith('#')) return;

  try {
    const feature = { attributes: {} };

    const columns = line.split(/\t/);

    feature['line'] = ++lineNo;
    feature['sequenceId'] = columns[0];
    feature['sequenceNumber'] = sequences[columns[0]].sequenceNumber;
    feature['type'] = columns[2];
    feature['startPosition'] = parseInt(columns[3]);
    feature['endPosition'] = parseInt(columns[4]);

    const infos = columns[8].split(';');
    if (infos && infos.length > 0) {
      const ids = infos.filter(id => id.startsWith('ID='));
      if (ids && ids.length > 0) {
        feature['id'] = ids[0].split('=')[1];
      }

      infos
        .filter(i => 'A' <= i[0][0] && i[0][0] <= 'Z')
        .forEach(i => {
          const is = i.split('=');
          feature['attributes'][is[0]] = is[1];
        });
    }

    console.log(feature);
  } catch (err) {
    console.error(line);
    console.error(err);
    process.exit(1);
  }
}

function processFeatures() {
  var buffer = '';
  var rs = fs.createReadStream(GFF);

  rs.on('data', function (chunk) {
    var lines = (buffer + chunk).split(/\r?\n/g);
    buffer = lines.pop();
    for (var i = 0; i < lines.length; ++i) {
      processFeature('' + lines[i]);
    }
  });

  rs.on('end', function () {
    console.log(`Finished - processed ${lineNo} segments`);
  });
}

function processSequence(line) {
  const sequenceRegex = />([^ ]*) (.*)/;
  const chromosomeRegex = /Homo sapiens chromosome ([0-9XY]+), GRCh38.p14 Primary Assembly/;

  const [, seqid, description] = sequenceRegex.exec(line);
  let chromosomeMatch = chromosomeRegex.exec(description);

  const sequence =
    chromosomeMatch !== null
      ? {
          sequenceNumber,
          seqid,
          isPrimary: true,
          chromosome: 'chr' + chromosomeMatch[1],
          description
        }
      : {
          sequenceNumber,
          seqid,
          isPrimary: false,
          description
        };

  sequenceNumber++;
  sequences[seqid] = sequence;
}

function processSequences() {
  var buffer = '';
  var rs = fs.createReadStream(FASTA);

  rs.on('data', function (chunk) {
    var lines = (buffer + chunk).split(/\r?\n/g);
    buffer = lines.pop();
    for (var i = 0; i < lines.length; ++i) {
      if (lines[i].startsWith('>')) processSequence(lines[i]);
    }
  });

  rs.on('end', function () {
    console.log(`Finished - processed ${sequenceNumber - 1024} sequences`);
    console.log(sequences);
  });
}

processSequences();
processFeatures();
