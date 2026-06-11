#!/usr/bin/env bash
set -u -o pipefail

# HTT LoMA local assembly driver/fallback.
# It prepares LoMA-compatible input directories, attempts to find/install LoMA,
# runs LoMA on the combined local reads and optional haplotype-split reads, and
# records enough diagnostic information to explain blockers in restricted environments.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${1:-${ROOT_DIR}/loma_output}"
INPUT_DIR="${ROOT_DIR}"
LOMA_REPO="https://github.com/kolikem/loma.git"
LOMA_TAG="v1.1.3"
LOMA_INSTALL_DIR="${ROOT_DIR}/tools/loma"
LOG_FILE="${OUT_DIR}/run.log"
SUMMARY_FILE="${OUT_DIR}/RUN_SUMMARY.md"

mkdir -p "${OUT_DIR}"
: > "${LOG_FILE}"

log() {
  printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" | tee -a "${LOG_FILE}"
}

record_summary_header() {
  cat > "${SUMMARY_FILE}" <<SUMMARY
# HTT LoMA run summary

- Target region: \`chr4:3074400-3075400\`
- LoMA package/tool: \`kolikem/loma\` (Localized Merging and Assembly), release \`${LOMA_TAG}\` when installable.
- LoMA input format: an **absolute input directory** containing one or more \`*.fastq\` files.
- Output root: \`${OUT_DIR}\`
- Log: \`run.log\`

SUMMARY
}

append_summary() {
  printf '%s\n' "$*" >> "${SUMMARY_FILE}"
}

copy_if_present() {
  local src="$1"
  local dest_dir="$2"
  mkdir -p "${dest_dir}"
  if [[ -s "${src}" ]]; then
    cp -f "${src}" "${dest_dir}/"
    log "Prepared input: ${src} -> ${dest_dir}/"
    return 0
  fi
  log "Missing or empty input: ${src}"
  return 1
}

find_loma() {
  if [[ -n "${LOMA_BIN:-}" && -x "${LOMA_BIN}" ]]; then
    printf '%s\n' "${LOMA_BIN}"
    return 0
  fi
  if command -v loma >/dev/null 2>&1; then
    command -v loma
    return 0
  fi
  if [[ -x "${LOMA_INSTALL_DIR}/loma" ]]; then
    printf '%s\n' "${LOMA_INSTALL_DIR}/loma"
    return 0
  fi
  return 1
}

attempt_loma_install() {
  mkdir -p "$(dirname "${LOMA_INSTALL_DIR}")"
  if [[ -d "${LOMA_INSTALL_DIR}/.git" ]]; then
    log "Existing LoMA checkout found at ${LOMA_INSTALL_DIR}"
    return 0
  fi
  log "Attempting local LoMA install: git clone --branch ${LOMA_TAG} ${LOMA_REPO} ${LOMA_INSTALL_DIR}"
  if git clone --depth 1 --branch "${LOMA_TAG}" "${LOMA_REPO}" "${LOMA_INSTALL_DIR}" >>"${LOG_FILE}" 2>&1; then
    chmod +x "${LOMA_INSTALL_DIR}/loma" || true
    log "LoMA cloned successfully."
    return 0
  fi
  log "LoMA clone failed; see run.log."
  return 1
}

check_dependencies() {
  local missing=0
  for exe in minimap2 mafft python3; do
    if command -v "${exe}" >/dev/null 2>&1; then
      log "Dependency found: ${exe} ($(command -v "${exe}"))"
    else
      log "Dependency missing: ${exe}"
      missing=1
    fi
  done
  return "${missing}"
}

run_loma_case() {
  local label="$1"
  local fastq_name="$2"
  local run_input="${OUT_DIR}/inputs/${label}"
  local run_output="${OUT_DIR}/${label}"
  local source_fastq="${INPUT_DIR}/${fastq_name}"

  mkdir -p "${run_input}" "${run_output}"
  if ! copy_if_present "${source_fastq}" "${run_input}"; then
    append_summary "## ${label}"
    append_summary ""
    append_summary "Status: **not run**; required input \`${fastq_name}\` was not present in \`${INPUT_DIR}\`."
    append_summary ""
    return 2
  fi

  if [[ -z "${LOMA_CMD:-}" ]]; then
    append_summary "## ${label}"
    append_summary ""
    append_summary "Status: **not run**; LoMA executable was unavailable. Prepared input directory: \`inputs/${label}\`."
    append_summary ""
    return 3
  fi

  log "Running LoMA for ${label}: ${LOMA_CMD} -I ${run_input} -O ${run_output}"
  append_summary "## ${label}"
  append_summary ""
  append_summary '```bash'
  append_summary "${LOMA_CMD} -I ${run_input} -O ${run_output}"
  append_summary '```'

  if bash "${LOMA_CMD}" -I "${run_input}" -O "${run_output}" >>"${LOG_FILE}" 2>&1; then
    log "LoMA completed for ${label}."
    append_summary ""
    append_summary "Status: **succeeded**."
  else
    log "LoMA failed for ${label}; see run.log."
    append_summary ""
    append_summary "Status: **failed**; see \`run.log\`."
    return 4
  fi

  local consensus_count=0
  if [[ -d "${run_output}/CONSENSUS" ]]; then
    consensus_count=$(find "${run_output}/CONSENSUS" -maxdepth 1 -type f -name '*.cs' | wc -l | tr -d ' ')
  fi
  append_summary "Consensus files generated: ${consensus_count}."

  if [[ "${label}" == "local_reads" ]]; then
    mapfile -t cs_files < <(find "${run_output}/CONSENSUS" -maxdepth 1 -type f -name '*.cs' | sort)
    if [[ "${#cs_files[@]}" -eq 1 ]]; then
      cp -f "${cs_files[0]}" "${OUT_DIR}/HTT.consensus.fasta"
      append_summary "One consensus sequence was produced; haplotype separation was **not achieved** for the combined-read run."
      append_summary "Saved: \`HTT.consensus.fasta\`."
    elif [[ "${#cs_files[@]}" -ge 2 ]]; then
      cp -f "${cs_files[0]}" "${OUT_DIR}/HTT.candidate_hap1.fasta"
      cp -f "${cs_files[1]}" "${OUT_DIR}/HTT.candidate_hap2.fasta"
      append_summary "Two or more candidate consensus sequences were produced."
      append_summary "Saved first two candidates as \`HTT.candidate_hap1.fasta\` and \`HTT.candidate_hap2.fasta\`."
    else
      append_summary "No consensus sequence was produced."
    fi
  fi
  append_summary ""
}

record_summary_header
log "Starting HTT LoMA preparation/run."

if ! find_loma >/tmp/htt_loma_cmd.txt; then
  attempt_loma_install || true
fi
LOMA_CMD="$(find_loma 2>/dev/null || true)"
if [[ -n "${LOMA_CMD}" ]]; then
  log "LoMA executable: ${LOMA_CMD}"
else
  log "LoMA executable unavailable."
fi

if ! check_dependencies; then
  log "One or more dependencies are missing. LoMA can only run after minimap2, MAFFT, and python3 are available."
fi

run_loma_case "local_reads" "HTT.local_reads_from_bam.fastq" || true
run_loma_case "hap1_reads" "HTT.hap1.reads.fastq" || true
run_loma_case "hap2_reads" "HTT.hap2.reads.fastq" || true

append_summary "## Generated files"
append_summary ""
find "${OUT_DIR}" -maxdepth 3 -type f | sort | sed "s#${OUT_DIR}/#- #" >> "${SUMMARY_FILE}"

log "Finished. Summary written to ${SUMMARY_FILE}."
