import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

const SCRIPT_DIR = path.dirname(fileURLToPath(import.meta.url));
const PACKAGE_ROOT = path.resolve(SCRIPT_DIR, "..", "..", "..");
const EXPERIMENT_DIR = path.join(PACKAGE_ROOT, "experiments", "08_pseudo_bulk_rnaseq_benchmark");

function parseArgs(argv) {
  const args = {
    inputDir: path.join(EXPERIMENT_DIR, "generated", "final_two_table_audit"),
    outDir: path.join(EXPERIMENT_DIR, "generated", "report"),
  };
  for (let i = 0; i < argv.length; i++) {
    const key = argv[i];
    const value = argv[i + 1];
    if (key === "--input-dir") {
      args.inputDir = value;
      i++;
    } else if (key === "--out-dir") {
      args.outDir = value;
      i++;
    } else if (key === "-h" || key === "--help") {
      console.log("Usage: node summarize_figure4_evidence_s18_s19.mjs [--input-dir DIR] [--out-dir DIR]");
      process.exit(0);
    } else {
      throw new Error(`Unknown argument: ${key}`);
    }
  }
  return args;
}

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let inQuotes = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    const next = text[i + 1];
    if (inQuotes) {
      if (ch === '"' && next === '"') {
        field += '"';
        i++;
      } else if (ch === '"') {
        inQuotes = false;
      } else {
        field += ch;
      }
    } else if (ch === '"') {
      inQuotes = true;
    } else if (ch === ",") {
      row.push(field);
      field = "";
    } else if (ch === "\n") {
      row.push(field);
      rows.push(row);
      row = [];
      field = "";
    } else if (ch !== "\r") {
      field += ch;
    }
  }
  if (field.length || row.length) {
    row.push(field);
    rows.push(row);
  }
  const [rawHeaders, ...data] = rows.filter((r) => r.some((v) => v !== ""));
  const headers = rawHeaders.map((h) => String(h).replace(/^\uFEFF/, "").trim());
  return data.map((r) => Object.fromEntries(headers.map((h, i) => [h, r[i] ?? ""])));
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

function toInt(value) {
  const n = Number.parseInt(value, 10);
  return Number.isFinite(n) ? n : 0;
}

function joinChanges(...parts) {
  return parts
    .flatMap((part) => String(part || "").split(";").map((x) => x.trim()))
    .filter((x) => x && x.toLowerCase() !== "none")
    .join("; ") || "None";
}

function s18Rows(sourceRows) {
  const strata = [
    "All defining sites observed",
    "Partially observed defining sites",
    "No defining site observed",
  ];
  const methods = [
    ["MiXCR-default", "MiXCR-default Exact Recovered"],
    ["FindAlleles", "FindAlleles Exact Recovered"],
    ["PanTCR", "PanTCR Exact Recovered"],
  ];
  const grouped = new Map();
  for (const row of sourceRows) {
    const key = `${row.Benchmark}\t${row.Dataset}`;
    if (!grouped.has(key)) grouped.set(key, []);
    grouped.get(key).push(row);
  }
  const out = [];
  for (const [key, rows] of grouped.entries()) {
    const [benchmark, dataset] = key.split("\t");
    const byStratum = new Map(rows.map((r) => [r["Observed Defining-site Stratum"], r]));
    const truthTotal = strata.reduce(
      (acc, s) => acc + toInt(byStratum.get(s)?.["No. of Non-default Truth Alleles"]),
      0,
    );
    for (const [method, column] of methods) {
      const recovered = {};
      for (const stratum of strata) {
        recovered[stratum] = toInt(byStratum.get(stratum)?.[column]);
      }
      const recoveredTotal = strata.reduce((acc, s) => acc + recovered[s], 0);
      out.push({
        Benchmark: benchmark,
        Dataset: dataset,
        Method: method,
        "NO. of Non-default Truth Alleles": truthTotal,
        "NO. of Recovered Alleles": recoveredTotal,
        "NO. of Not Recovered Alleles": truthTotal - recoveredTotal,
        "Recovered with All Changes Observed": recovered["All defining sites observed"],
        "Recovered with Partial Changes Observed": recovered["Partially observed defining sites"],
        "Recovered with No Change Observed": recovered["No defining site observed"],
      });
    }
  }
  return out;
}

function s19Rows(sourceRows) {
  return sourceRows.map((row) => ({
    Dataset: row.Dataset,
    Gene: row.Gene,
    "Truth Allele": row["Truth allele"],
    "Default-relative Change": row["Change Relative to Default Allele"] || "None",
    "Evidence-supported Change": joinChanges(
      row["Evidence-supported Defining Site"],
      row["Mixed Observed Defining Site"],
    ),
    "Graph-imputed Change": joinChanges(row["Graph-imputed Defining Site"]),
    "Inferred Sequence": row["Inferred Sequence"],
  }));
}

const { inputDir, outDir } = parseArgs(process.argv.slice(2));
await fs.mkdir(outDir, { recursive: true });

const methodRows = parseCsv(
  await fs.readFile(path.join(inputDir, "method_comparison_by_observed_defining_site.csv"), "utf8"),
);
const exampleRows = parseCsv(
  await fs.readFile(path.join(inputDir, "pantcr_recovered_scbulk_examples.csv"), "utf8"),
);

const s18 = s18Rows(methodRows);
const s19 = s19Rows(exampleRows);
const s18Out = path.join(outDir, "supplementary_table_s18_source_rows.csv");
const s19Out = path.join(outDir, "supplementary_table_s19_source_rows.csv");
const noteOut = path.join(outDir, "figure4e_s18_s19_notes.md");

await fs.writeFile(
  s18Out,
  writeCsv(s18, [
    "Benchmark",
    "Dataset",
    "Method",
    "NO. of Non-default Truth Alleles",
    "NO. of Recovered Alleles",
    "NO. of Not Recovered Alleles",
    "Recovered with All Changes Observed",
    "Recovered with Partial Changes Observed",
    "Recovered with No Change Observed",
  ]),
);
await fs.writeFile(
  s19Out,
  writeCsv(s19, [
    "Dataset",
    "Gene",
    "Truth Allele",
    "Default-relative Change",
    "Evidence-supported Change",
    "Graph-imputed Change",
    "Inferred Sequence",
  ]),
);
await fs.writeFile(
  noteOut,
  [
    "# Figure 4E / Supplementary Tables S18-S19 Source Summary",
    "",
    "S18 uses unfiltered MiXCR-derived observed-region coverage for method-comparison strata.",
    "S19 uses retained PanTCR evidence to separate evidence-supported and graph-imputed recovered changes.",
    "These files are compact source summaries and do not recreate the formal supplementary workbook layout.",
    "",
  ].join("\n"),
);

console.log(s18Out);
console.log(s19Out);
console.log(noteOut);
