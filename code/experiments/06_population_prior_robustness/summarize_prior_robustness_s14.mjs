import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

const SCRIPT_DIR = path.dirname(fileURLToPath(import.meta.url));
const PACKAGE_ROOT = path.resolve(SCRIPT_DIR, "..", "..", "..");
const EXPERIMENT_DIR = path.join(PACKAGE_ROOT, "experiments", "06_population_prior_robustness");

function parseArgs(argv) {
  const args = {
    sourceCsv: path.join(EXPERIMENT_DIR, "generated", "summary", "population_prior_robustness_s14_source.csv"),
    outDir: path.join(EXPERIMENT_DIR, "generated", "report"),
  };
  for (let i = 0; i < argv.length; i++) {
    const key = argv[i];
    const value = argv[i + 1];
    if (key === "--source-csv") {
      args.sourceCsv = value;
      i++;
    } else if (key === "--out-dir") {
      args.outDir = value;
      i++;
    } else if (key === "-h" || key === "--help") {
      console.log("Usage: node summarize_prior_robustness_s14.mjs [--source-csv CSV] [--out-dir DIR]");
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

const { sourceCsv, outDir } = parseArgs(process.argv.slice(2));
await fs.mkdir(outDir, { recursive: true });
const records = parseCsv(await fs.readFile(sourceCsv, "utf8"));

const columns = [
  "Setting",
  "Mode",
  "NO. of Truth Alleles",
  "NO. of Predicted Alleles",
  "TP",
  "FP",
  "FN",
  "Precision",
  "Recall",
  "F1-score",
];
const outCsv = path.join(outDir, "population_prior_robustness_s14_table.csv");
const outNote = path.join(outDir, "population_prior_robustness_s14_notes.md");

await fs.writeFile(outCsv, writeCsv(records, columns));
await fs.writeFile(
  outNote,
  [
    "# Population-Prior Robustness",
    "",
    "Scope: AFR TRBV allele inference in Scenario A and Scenario B in silico benchmarks.",
    "Modes: no graph prior, EUR-mismatched prior, global prior, and AFR-matched prior.",
    "Report files contain compact source values for the population-prior robustness result summary.",
    "",
  ].join("\n"),
);

console.log(outCsv);
console.log(outNote);
