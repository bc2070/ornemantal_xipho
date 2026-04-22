#!/bin/bash
set -e

WORK_DIR=$1
MARKER_FILENAME=$2
IND_NUM=$3
S_VALUE=$4

if [ -z "$WORK_DIR" ] || [ -z "$MARKER_FILENAME" ] || [ -z "$IND_NUM" ] || [ -z "$S_VALUE" ]; then
    echo "Usage: $0 <working_directory> <marker_filename> <individual_number> <S_value>" >&2
    exit 1
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SUMMARY_FILE="${SCRIPT_DIR}/Ss.hm.s.txt"

SC1=/public3/group_crf/home/b20gengjlin3/Xiphophorus/getfind.py
SC2=/fast3/group_crf/home/b20gengjlin3/xipho_simul/neutral_seg/APS/tmp/01.batch.product.4.2.py
SC3=/fast3/group_crf/home/b20gengjlin3/xipho_simul/neutral_seg/APS/tmp/extra_analysis.py
CHROM_FILE=/public3/group_crf/home/b20gengjlin3/Xiphophorus/chrom.txt
H_VALUE_FILE=/fast3/group_crf/home/b20gengjlin3/xipho_simul/neutral_seg/h_value_tsv/hm.h_value.tsv
MARKER_POS_FILE=/fast3/group_crf/home/b20gengjlin3/xipho_simul/neutral_seg/ahmm_marker/APS1.marker.pos.txt

cd "$WORK_DIR"

MARKER_FULL_PATH="./$MARKER_FILENAME"
BASENAME=${MARKER_FILENAME%_markers.txt}
LINE_COUNT=$(wc -l < "$MARKER_FULL_PATH")

python3 "$SC1" "$MARKER_FULL_PATH" "$BASENAME" "$IND_NUM" "$LINE_COUNT" "$MARKER_POS_FILE" > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "FATAL: getfind.py failed." >&2
    exit 1
fi

wait_timeout=600; wait_counter=0; expected_files=$IND_NUM; stability_sleep=3
while true; do
    if [[ $wait_counter -ge $wait_timeout ]]; then
        echo "Error: Timed out after ${wait_timeout}s waiting for find files." >&2; exit 1
    fi
    find_files=( "${BASENAME}"_find*.txt )
    if [[ ! -f "${find_files[0]}" ]]; then file_count=0; else file_count=${#find_files[@]}; fi
    if [[ $file_count -ge $expected_files ]]; then
        current_size=$(stat -c %s "${BASENAME}"_find*.txt | awk '{total+=$1} END{print total}')
        sleep $stability_sleep
        new_size=$(stat -c %s "${BASENAME}"_find*.txt | awk '{total+=$1} END{print total}')
        if [[ "$current_size" == "$new_size" ]] && [[ -n "$current_size" ]]; then
            break
        fi
    fi
    sleep 1; ((wait_counter++))
done

find_files=( "${BASENAME}"_find*.txt )
OUTPUT_BED="${BASENAME}.lk.bed"
python3 "$SC2" "$CHROM_FILE" "$H_VALUE_FILE" "${find_files[@]}" > "$OUTPUT_BED"

if [ ! -s "$OUTPUT_BED" ]; then
    echo "Error: Output file '$OUTPUT_BED' is empty or was not created." >&2
    exit 1
fi

ANALYSIS_RESULTS=$(python3 "$SC3" "$OUTPUT_BED" "$BASENAME" "$CHROM_FILE")
if [ $? -ne 0 ]; then
    echo "FATAL: extra_analysis.py failed." >&2
    exit 1
fi

read -r mean_len_no14 stdev_len_no14 mean_proportion eff_len <<< "$ANALYSIS_RESULTS"

IFS='_' read -r p1 p2 p3 p4 <<< "$BASENAME"

OUTPUT_LINE="${p1}\t${p2}\t${p3}\t${p4}\t${S_VALUE}\t${mean_len_no14}\t${stdev_len_no14}\t${mean_proportion}\t${eff_len}"

echo -e "$OUTPUT_LINE" >> "$SUMMARY_FILE"

rm -f "${BASENAME}"_find*.txt
rm -f "$MARKER_FULL_PATH"
