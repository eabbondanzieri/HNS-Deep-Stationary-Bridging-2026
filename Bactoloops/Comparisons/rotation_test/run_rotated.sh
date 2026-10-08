#!/bin/bash
# Run every caller (settings as run_callers.sh) and the ImageJ macro on one rotated input.
# Usage: bash rotation_test/run_rotated.sh <dataset>_r<k>
set -u
cd "$(dirname "$0")"
HX="${HICEXPLORER_BIN:-$(dirname "$(command -v hicDetectLoops)")}"
LC="${LOOPCOMPARE_BIN:-$(dirname "$(command -v cooltools)")}"
PK="${PEAKACHU_BIN:-$(dirname "$(command -v peakachu)")}"
T=../tools; JT=$T/juicer_tools_1.22.01.jar
d=$1; in=inputs/$d.cool; len=$(cut -f2 inputs/$d.chrom.sizes); chr=$(cut -f1 inputs/$d.chrom.sizes)
o=runs; mkdir -p $o/{cooltools,chromosight,mustache,peakachu,hiccups,sip,hic,logs}
"${FIJI:-fiji}" --allow-multiple --no-splash -- --headless -macro "$(pwd)/macros/$d.ijm" > $o/logs/$d.bactoloops.log 2>&1
$HX/cooltools expected-cis $in -o $o/cooltools/$d.expected.tsv > $o/logs/$d.cooltools.log 2>&1
$HX/cooltools dots $in $o/cooltools/$d.expected.tsv --max-loci-separation $len -o $o/cooltools/$d.tsv >> $o/logs/$d.cooltools.log 2>&1
$HX/chromosight detect --pattern loops --max-dist $len --threads 4 $in $o/chromosight/$d > $o/logs/$d.chromosight.log 2>&1
$LC/mustache -f $in -r 10000 -d $len -p 4 -o $o/mustache/$d.tsv > $o/logs/$d.mustache.log 2>&1
$PK/peakachu score_genome -r 10000 -p $in -C $chr --clr-weight-name weight -O $o/peakachu/$d.scores.bedpe -m $T/peakachu_models/high-confidence.2billion.10kb.w6.pkl -u $(( len / 10000 + 1 )) > $o/logs/$d.peakachu.log 2>&1
$PK/peakachu pool -r 10000 -i $o/peakachu/$d.scores.bedpe -o $o/peakachu/$d.tsv -t 0.9 >> $o/logs/$d.peakachu.log 2>&1
awk '{print $1"\t"$2+10000000}' inputs/$d.chrom.sizes > $o/hic/$d.padded.sizes
rm -f $o/hic/$d.hic; java -Xmx4g -jar $JT pre -r 10000 inputs/$d.short.txt.gz $o/hic/$d.hic $o/hic/$d.padded.sizes > $o/logs/$d.pre.log 2>&1
rm -rf $o/hiccups/$d; java -Xmx4g -jar $JT hiccups --cpu --ignore-sparsity -r 10000 $o/hic/$d.hic $o/hiccups/$d > $o/logs/$d.hiccups.log 2>&1
rm -rf $o/sip/$d; java -Djava.awt.headless=true -Xmx4g -jar $T/SIP-1.6.5.jar hic $o/hic/$d.hic $o/hic/$d.padded.sizes $o/sip/$d $T/juicer_tools_1.13.02.jar -res 10000 > $o/logs/$d.sip.log 2>&1
echo "$d done"
