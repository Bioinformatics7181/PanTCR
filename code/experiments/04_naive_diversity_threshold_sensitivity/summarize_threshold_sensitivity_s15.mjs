import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

const SCRIPT_DIR = path.dirname(fileURLToPath(import.meta.url));
const PACKAGE_ROOT = path.resolve(SCRIPT_DIR, "..", "..", "..");
const EXPERIMENT_DIR = path.join(PACKAGE_ROOT, "experiments", "04_naive_diversity_threshold_sensitivity");

function parseArgs(argv) {
  const args = {
    summaryDir: path.join(EXPERIMENT_DIR, "generated", "summary"),
    outputDir: path.join(EXPERIMENT_DIR, "generated", "report"),
  };
  for (let i = 0; i < argv.length; i++) {
    const key = argv[i];
    const value = argv[i + 1];
    if (key === "--summary-dir") {
      args.summaryDir = value;
      i++;
    } else if (key === "--out-dir") {
      args.outputDir = value;
      i++;
    } else if (key === "-h" || key === "--help") {
      console.log("Usage: node summarize_threshold_sensitivity_s15.mjs [--summary-dir DIR] [--out-dir DIR]");
      process.exit(0);
    } else {
      throw new Error(`Unknown argument: ${key}`);
    }
  }
  return args;
}

function parseCsvLine(line) {
  const out = [];
  let cur = "";
  let inQuote = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (ch === '"') {
      if (inQuote && line[i + 1] === '"') {
        cur += '"';
        i++;
      } else {
        inQuote = !inQuote;
      }
    } else if (ch === "," && !inQuote) {
      out.push(cur);
      cur = "";
    } else {
      cur += ch;
    }
  }
  out.push(cur);
  return out;
}

function parseCsv(text) {
  const lines = text.trim().split(/\r?\n/).filter(Boolean);
  const header = parseCsvLine(lines.shift() || "");
  return lines.map((line) => {
    const cells = parseCsvLine(line);
    return Object.fromEntries(header.map((h, i) => [h.replace(/^\uFEFF/, ""), cells[i] ?? ""]));
  });
}

function csvEscape(value) {
  const text = String(value ?? "");
  return /[",\r\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
}

function writeCsv(rows, columns) {
  return [
    columns.join(","),
    ...rows.map((row) => columns.map((col) => csvEscape(row[col])).join(",")),
  ].join("\n") + "\n";
}

function toDisplayMethod(method) {
  if (method === "PanTCR_fixed_graph") return "PanTCR";
  if (method === "PanTCR_NP_no_prior") return "PanTCR-NP";
  return method;
}

const { summaryDir, outputDir } = parseArgs(process.argv.slice(2));
await fs.mkdir(outputDir, { recursive: true });

const metricsCsv = await fs.readFile(
  path.join(summaryDir, "threshold_precision_recall_f1_by_gene.csv"),
  "utf8",
);
const deltaCsv = await fs.readFile(
  path.join(summaryDir, "threshold_f1_delta_summary.csv"),
  "utf8",
);

const metricRows = parseCsv(metricsCsv)
  .map((row) => ({
    Gene: row.gene,
    Method: toDisplayMethod(row.method_label),
    "Inference min_naive": Number(row.inference_min_naive),
    Precision: row.Precision,
    Recall: row.Recall,
    "F1-score": row["F1-score"],
  }))
  .sort((a, b) => a.Gene.localeCompare(b.Gene) || a["Inference min_naive"] - b["Inference min_naive"] || a.Method.localeCompare(b.Method));

const deltaRows = parseCsv(deltaCsv)
  .map((row) => ({
    Gene: row.gene,
    "Inference min_naive": Number(row.inference_min_naive),
    "PanTCR F1": row.PanTCR_fixed_graph,
    "PanTCR-NP F1": row.PanTCR_NP_no_prior,
    "PanTCR minus PanTCR-NP F1": row.PanTCR_minus_NP_F1,
  }))
  .sort((a, b) => a.Gene.localeCompare(b.Gene) || a["Inference min_naive"] - b["Inference min_naive"]);

const metricOut = path.join(outputDir, "threshold_sensitivity_main_results.csv");
const deltaOut = path.join(outputDir, "threshold_sensitivity_f1_delta.csv");
const noteOut = path.join(outputDir, "threshold_sensitivity_notes.md");

await fs.writeFile(metricOut, writeCsv(metricRows, ["Gene", "Method", "Inference min_naive", "Precision", "Recall", "F1-score"]));
await fs.writeFile(deltaOut, writeCsv(deltaRows, ["Gene", "Inference min_naive", "PanTCR F1", "PanTCR-NP F1", "PanTCR minus PanTCR-NP F1"]));
await fs.writeFile(
  noteOut,
  [
    "# Inference-Only Diversity-Threshold Sensitivity",
    "",
    "Scope: Scenario A fixed graph; only the sample-level inference min_naive threshold is changed.",
    "Compared methods: PanTCR and PanTCR-NP under inference min_naive = 0, 1, 2, 3.",
    "Report files contain compact source values for the threshold-sensitivity result summary.",
    "",
  ].join("\n"),
);

console.log(metricOut);
console.log(deltaOut);
console.log(noteOut);
