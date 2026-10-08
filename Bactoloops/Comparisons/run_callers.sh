#!/bin/bash
# Run each published loop caller on the comparison inputs (build_inputs.py).
# Usage: bash run_callers.sh <tool> [dataset ...]      tools: cooltools chromosight mustache peakachu hiccups sip
# Defaults are each tool's own, except the maximum loop distance, which is opened up to the
# whole chromosome so every tool searches the same range of distances as the ImageJ macro
# (bacterial loops span the whole chromosome; the tools' 2-3 Mb defaults are tuned for mammals).
# None of these tools treats the chromosome as circular, so loops spanning oriC (E. coli) or
# the map origin are out of reach for all of them.
set -euo pipefail
cd "$(dirname "$0")"
HX="${HICEXPLORER_BIN:-$(dirname "$(command -v hicDetectLoops)")}"
LC="${LOOPCOMPARE_BIN:-$(dirname "$(command -v cooltools)")}"
PK="${PEAKACHU_BIN:-$(dirname "$(command -v peakachu)")}"
JT=tools/juicer_tools_1.22.01.jar
tool=$1; shift
datasets=${*:-wt96 wtexp hns96 bsub}

for d in $datasets; do
  in=inputs/$d.cool
  len=$(cut -f2 inputs/$d.chrom.sizes)
  out=runs/$tool; mkdir -p $out
  case $tool in
    cooltools)   # HiCCUPS-style donut/lowleft/horizontal/vertical kernels, FDR 0.1
      $HX/cooltools expected-cis $in -o $out/$d.expected.tsv
      $HX/cooltools dots $in $out/$d.expected.tsv --max-loci-separation $len -o $out/$d.tsv ;;
    chromosight) # template matching with the default loop kernel
      $HX/chromosight detect --pattern loops --max-dist $len --threads 4 $in $out/$d > $out/$d.log 2>&1 ;;
    mustache)    # scale-space blob detection; default thresholds (FDR 0.2)
      $LC/mustache -f $in -r 10000 -d $len -p 4 -o $out/$d.tsv > $out/$d.log 2>&1 || true ;;
    peakachu)    # random forest, pretrained high-confidence model, as `peakachu depth` suggests
      # `peakachu depth` scales depth to a human genome; every map here maps to its top model
      model=tools/peakachu_models/high-confidence.2billion.10kb.w6.pkl
      nb=$(( len / 10000 + 1 ))
      $PK/peakachu score_genome -r 10000 -p $in -C $(cut -f1 inputs/$d.chrom.sizes) --clr-weight-name weight -O $out/$d.scores.bedpe -m $model -u $nb > $out/$d.log 2>&1
      $PK/peakachu pool -r 10000 -i $out/$d.scores.bedpe -o $out/$d.tsv -t 0.9 >> $out/$d.log 2>&1 ;;
    hiccups)     # juicer HiCCUPS on the CPU, 10 kb defaults
      # The CPU version indexes one bin past the end of the map (ArrayIndexOutOfBounds at the
      # last bin) and drops that tile. The .hic is therefore built with 10 Mb of empty sequence
      # appended to the chromosome, so the failing tile holds no data (2 Mb and 10 Mb pads give
      # identical calls on WT 96 h). The padding carries no contacts, so no loop can be called in it.
      mkdir -p inputs/hic
      awk '{print $1"\t"$2+10000000}' inputs/$d.chrom.sizes > inputs/hic/$d.padded.sizes
      [ -f inputs/hic/$d.hic ] || java -Xmx4g -jar $JT pre -r 10000 inputs/$d.short.txt.gz inputs/hic/$d.hic inputs/hic/$d.padded.sizes > inputs/hic/$d.pre.log 2>&1
      java -Xmx4g -jar $JT hiccups --cpu --ignore-sparsity -r 10000 inputs/hic/$d.hic $out/$d > $out/$d.log 2>&1 || true ;;
    sip)         # SIP on the .hic built for HiCCUPS, default parameters. SIP reads the .hic through
                 # juicer_tools 1.13, as its authors recommend: 1.22 dumps only the first expected value
      java -Djava.awt.headless=true -Xmx4g -jar tools/SIP-1.6.5.jar hic inputs/hic/$d.hic inputs/hic/$d.padded.sizes $out/$d tools/juicer_tools_1.13.02.jar -res 10000 > $out/$d.log 2>&1 || true ;;
  esac
  echo "$tool $d done"
done
