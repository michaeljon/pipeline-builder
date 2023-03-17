const fs = require('fs');
const path = require('path');

if (process.argv.length < 4) {
  console.log('Usage: node ' + process.argv[1] + ' <BUCKET> <FILENAME>');
  process.exit(1);
}

var bucket = process.argv[2];
var file = process.argv[3];

var base = path.basename(file);
var timestamp = base.split('_')[0] + ' ' + base.split('_')[1];
var location = base.split('_')[2];
var sequence = location.replace('.json', '');

try {
  var row = [];
  var json = fs.readFileSync(file, 'utf8');

  var stats = JSON.parse(json.replace(/\bNaN\b/g, '"***NaN***"'), function (key, value) {
    return value === '***NaN***' ? NaN : value;
  });

  const sample = Object.keys(stats.report_data_sources.Samtools.stats)[0];

  row.push(sequence);
  row.push(sample);
  row.push(timestamp);
  row.push(bucket + '/' + location);

  // parse the json
  row.push((stats.report_saved_raw_data.multiqc_general_stats[sample + '.duplication_metrics']['biobambam2_mqc-generalstats-biobambam2-PERCENT_DUPLICATION'] * 100.0).toFixed(4));

  row.push(stats.report_saved_raw_data.multiqc_general_stats[sample]['Samtools_mqc-generalstats-samtools-flagstat_total']);
  row.push(stats.report_saved_raw_data.multiqc_general_stats[sample]['Samtools_mqc-generalstats-samtools-mapped_passed']);

  // special case, picard is reporting its output differently sometimes
  row.push(
    (
      (stats.report_saved_raw_data.multiqc_general_stats[sample]['Picard_mqc-generalstats-picard-PCT_PF_READS_ALIGNED'] ||
        stats.report_saved_raw_data.multiqc_general_stats[sample + '_0']['Picard_mqc-generalstats-picard-PCT_PF_READS_ALIGNED']) * 100.0
    ).toFixed(4)
  );
  row.push(
    (
      stats.report_saved_raw_data.multiqc_general_stats[sample]['Picard_mqc-generalstats-picard-MEDIAN_COVERAGE'] ||
      stats.report_saved_raw_data.multiqc_general_stats[sample + '_0']['Picard_mqc-generalstats-picard-MEDIAN_COVERAGE']
    ).toFixed(4)
  );
  row.push(
    (
      stats.report_saved_raw_data.multiqc_general_stats[sample]['Picard_mqc-generalstats-picard-MEAN_COVERAGE'] ||
      stats.report_saved_raw_data.multiqc_general_stats[sample + '_0']['Picard_mqc-generalstats-picard-MEAN_COVERAGE']
    ).toFixed(4)
  );
  row.push(
    (
      stats.report_saved_raw_data.multiqc_general_stats[sample]['Picard_mqc-generalstats-picard-SD_COVERAGE'] ||
      stats.report_saved_raw_data.multiqc_general_stats[sample + '_0']['Picard_mqc-generalstats-picard-SD_COVERAGE']
    ).toFixed(4)
  );
  row.push(
    (
      (stats.report_saved_raw_data.multiqc_general_stats[sample]['Picard_mqc-generalstats-picard-PCT_30X'] ||
        stats.report_saved_raw_data.multiqc_general_stats[sample + '_0']['Picard_mqc-generalstats-picard-PCT_30X']) * 100.0
    ).toFixed(4)
  );

  row.push(stats.report_saved_raw_data.multiqc_general_stats[sample]['Samtools_mqc-generalstats-samtools-error_rate'].toFixed(6));
  row.push(stats.report_saved_raw_data.multiqc_general_stats[sample]['Samtools_mqc-generalstats-samtools-non_primary_alignments'].toFixed(0));
  row.push(stats.report_saved_raw_data.multiqc_general_stats[sample]['Samtools_mqc-generalstats-samtools-reads_mapped'].toFixed(0));
  row.push(stats.report_saved_raw_data.multiqc_general_stats[sample]['Samtools_mqc-generalstats-samtools-reads_mapped_percent'].toFixed(4));
  row.push(stats.report_saved_raw_data.multiqc_general_stats[sample]['Samtools_mqc-generalstats-samtools-reads_properly_paired_percent'].toFixed(4));
  row.push(stats.report_saved_raw_data.multiqc_general_stats[sample]['Samtools_mqc-generalstats-samtools-reads_MQ0_percent'].toFixed(4));
  row.push(stats.report_saved_raw_data.multiqc_general_stats[sample]['Samtools_mqc-generalstats-samtools-raw_total_sequences'].toFixed(0));

  console.log(row.join('\t'));
} catch (e) {
  console.log('Processing: ' + file);
  console.log('Error:', e.stack);
  process.exit(1);
}