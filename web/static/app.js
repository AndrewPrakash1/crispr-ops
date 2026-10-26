const gene = document.querySelector("#gene");
const population = document.querySelector("#population");
const note = document.querySelector("#note");
const guideList = document.querySelector("#guide-list");
const trace = document.querySelector("#trace");
const lead = document.querySelector("#lead");
const qualifier = document.querySelector("#qualifier");
const verdicts = document.querySelector("#verdicts");
const table = document.querySelector("#table");
const siteCard = document.querySelector("#site");
const benchOut = document.querySelector("#bench-out");
const validOut = document.querySelector("#valid-out");
const limits = document.querySelector("#limits");
const runButton = document.querySelector("#run");
const savedRun = document.querySelector("#saved-run");
const loadRunButton = document.querySelector("#load-run");
const traceMode = document.querySelector("#trace-mode");
const gridStatus = document.querySelector("#grid-status");
const gridTable = document.querySelector("#grid-table");
const gridCsvButton = document.querySelector("#grid-csv");
const sitesCsvButton = document.querySelector("#sites-csv");
const gridJsonButton = document.querySelector("#grid-json");

const guides = [];
let panelSnapshot = null;

function field(value, placeholder) {
  const input = document.createElement("input");
  input.value = value;
  input.placeholder = placeholder;
  input.autocomplete = "off";
  return input;
}

function addGuide(id, spacer) {
  const row = document.createElement("div");
  row.className = "guide-row";
  const idInput = field(id, "id");
  const spacerInput = field(spacer, "20 nt spacer");
  row.append(idInput, spacerInput);
  guideList.append(row);
  guides.push({ idInput, spacerInput });
}

fetch("/api/guides").then((r) => r.json()).then((data) => {
  const panel = (data.guides && data.guides.HBB) || [];
  panel.forEach((guide) => addGuide(guide.id, guide.spacer));
  if (!panel.length) addGuide("guide-1", "");
});

fetch("/api/validation").then((r) => r.json()).then((data) => {
  validOut.textContent = data.available
    ? JSON.stringify(data, null, 2)
    : data.reason;
});

async function refreshSavedRuns() {
  const response = await fetch("/api/runs");
  if (!response.ok) throw new Error(`Could not list saved runs (${response.status}).`);
  const data = await response.json();
  savedRun.replaceChildren();
  const runs = data.runs || [];
  runs.forEach((run) => {
    const option = document.createElement("option");
    option.value = run.run_id;
    option.textContent = `${run.gene || "Unknown"} / ${run.population || "unknown"} · ${new Date(run.created_at).toLocaleString()}${run.finished ? "" : " · incomplete"}`;
    option.disabled = !run.finished;
    savedRun.append(option);
  });
  if (!savedRun.options.length) {
    const option = document.createElement("option");
    option.value = "";
    option.textContent = "No saved runs";
    savedRun.append(option);
  }
  loadRunButton.disabled = !runs.some((run) => run.finished);
}

refreshSavedRuns().catch((error) => {
  savedRun.replaceChildren();
  const option = document.createElement("option");
  option.value = "";
  option.textContent = "Saved runs unavailable";
  savedRun.append(option);
  loadRunButton.disabled = true;
  traceMode.textContent = String(error);
});

async function loadPanelSnapshot() {
  const response = await fetch("/api/panel");
  if (!response.ok) throw new Error(`Frozen panel snapshot unavailable (${response.status}).`);
  panelSnapshot = await response.json();
  const rows = panelSnapshot.grid || [];
  const sites = panelSnapshot.sites || [];
  const evoRows = sites.filter((row) => row.evo2_direction).length;
  gridStatus.textContent = `Frozen local snapshot · ${rows.length} guide/population rows · ${sites.length} site rows · ${evoRows} Evo 2 directional outputs · loading this snapshot makes no live calls.`;
  const columns = ["gene", "population", "guide_id", "reference_verdict", "population_verdict", "flipped", "escalated", "n_sites"];
  const head = document.createElement("tr");
  columns.forEach((column) => {
    const th = document.createElement("th");
    th.textContent = column.replaceAll("_", " ");
    head.append(th);
  });
  gridTable.replaceChildren(head);
  rows.forEach((row) => {
    const tr = document.createElement("tr");
    columns.forEach((column) => {
      const td = document.createElement("td");
      td.textContent = row[column] == null ? "" : String(row[column]);
      tr.append(td);
    });
    gridTable.append(tr);
  });
  gridCsvButton.disabled = !rows.length;
  sitesCsvButton.disabled = !sites.length;
  gridJsonButton.disabled = !rows.length;
}

loadPanelSnapshot().catch((error) => {
  gridStatus.textContent = String(error);
});

function csvText(rows) {
  if (!rows.length) return "";
  const columns = Object.keys(rows[0]);
  const quote = (value) => `"${String(value ?? "").replaceAll('"', '""')}"`;
  return [columns.map(quote).join(","), ...rows.map((row) => columns.map((column) => quote(row[column])).join(","))].join("\r\n");
}

function downloadFile(filename, content, type) {
  const url = URL.createObjectURL(new Blob([content], { type }));
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = filename;
  anchor.click();
  URL.revokeObjectURL(url);
}

gridCsvButton.addEventListener("click", () => {
  if (panelSnapshot) downloadFile("grid.csv", csvText(panelSnapshot.grid || []), "text/csv;charset=utf-8");
});

sitesCsvButton.addEventListener("click", () => {
  if (panelSnapshot) downloadFile("grid_sites.csv", csvText(panelSnapshot.sites || []), "text/csv;charset=utf-8");
});

gridJsonButton.addEventListener("click", () => {
  if (panelSnapshot) downloadFile("grid.json", JSON.stringify(panelSnapshot, null, 2), "application/json");
});

function pill(stage) {
  const span = document.createElement("span");
  span.className = `pill ${stage}`;
  span.textContent = stage;
  return span;
}

function addTrace(event) {
  const item = document.createElement("li");
  const title = document.createElement("div");
  title.append(pill(event.stage), document.createTextNode(event.title || ""));
  item.append(title);
  if (event.detail) {
    const detail = document.createElement("p");
    detail.className = "detail";
    detail.textContent = event.detail;
    item.append(detail);
  }
  trace.append(item);
  trace.scrollTop = trace.scrollHeight;
}

function bases(seq, other) {
  const row = document.createElement("div");
  row.className = "seq";
  [...(seq || "")].forEach((base, i) => {
    const span = document.createElement("span");
    span.textContent = base;
    span.className = `base${base}`;
    if (other && other[i] && other[i] !== base) span.classList.add("mm");
    row.append(span);
  });
  return row;
}

function showSite(row) {
  siteCard.replaceChildren();
  const title = document.createElement("p");
  title.textContent = `${row.guide_id}  ${row.chrom}:${row.pos}${row.strand}  ${row.variant || ""}  AF ${row.af ?? "—"}`;
  siteCard.append(title);
  siteCard.append(bases(row.ref_window, row.alt_window));
  siteCard.append(bases(row.alt_window, row.ref_window));
  const scores = document.createElement("p");
  scores.textContent = `CFD ${Number(row.cfd_ref).toFixed(4)} → ${Number(row.cfd_alt).toFixed(4)}   Δ ${Number(row.delta).toFixed(4)}   mm ${row.mm_ref} → ${row.mm_alt}`;
  siteCard.append(scores);
  const reading = document.createElement("p");
  reading.textContent = row.reading;
  siteCard.append(reading);
  const evo = document.createElement("p");
  const opinion = row.evo2 || {};
  if (opinion.direction) {
    evo.textContent = `Evo 2 and CFD move in the ${opinion.direction}. ${opinion.reading || ""} Δlogp ${Number(opinion.evo_delta).toExponential(2)}.`;
  } else {
    const reason = (opinion.ref && opinion.ref.reason) || opinion.reason || "Evo 2 was not called for this site.";
    evo.textContent = `Evo 2: ${reason}`;
  }
  siteCard.append(evo);
  if (row.pairs && row.pairs.length) {
    const pairs = document.createElement("p");
    pairs.textContent = row.pairs.map((pair) => (
      `pos ${pair.position} ${pair.ref_base}→${pair.alt_base ?? pair.ref_base} weight ${Number(pair.ref_weight).toFixed(3)}→${pair.alt_weight == null ? "—" : Number(pair.alt_weight).toFixed(3)}`
    )).join(" · ");
    siteCard.append(pairs);
  }
}

function renderResult(payload) {
  if (payload.refused) {
    lead.textContent = payload.refusal ? payload.refusal.reason : "Refused.";
    qualifier.textContent = payload.qualifier || "";
    return;
  }
  lead.textContent = payload.lead || "";
  qualifier.textContent = payload.qualifier || "";
  verdicts.replaceChildren();
  (payload.guides || []).forEach((guide) => {
    const box = document.createElement("div");
    box.className = "verdict";
    box.textContent = `${guide.id}: reference ${guide.reference_verdict}, population ${guide.population_verdict}${guide.flipped ? ", flipped" : ""}${guide.escalated ? ", PAM expanded" : ""}`;
    verdicts.append(box);
  });
  table.replaceChildren();
  const head = document.createElement("tr");
  ["Guide", "Site", "mm", "CFD ref", "CFD alt", "Δ", "AF", "Reading"].forEach((label) => {
    const cell = document.createElement("th");
    cell.textContent = label;
    head.append(cell);
  });
  table.append(head);
  const rows = payload.table || [];
  if (!rows.length) {
    const empty = document.createElement("tr");
    const cell = document.createElement("td");
    cell.colSpan = 8;
    cell.textContent = "No off-target window in this run carries a common variant that changes CFD.";
    empty.append(cell);
    table.append(empty);
  }
  rows.forEach((row) => {
    const tr = document.createElement("tr");
    tr.className = "pick";
    const cells = [
      row.guide_id,
      `${row.chrom}:${row.pos}${row.strand}`,
      `${row.mm_ref}→${row.mm_alt}`,
      Number(row.cfd_ref).toFixed(4),
      Number(row.cfd_alt).toFixed(4),
      Number(row.delta).toFixed(4),
      row.af == null ? "—" : Number(row.af).toFixed(3),
      row.reading,
    ];
    cells.forEach((value) => {
      const td = document.createElement("td");
      td.textContent = value;
      tr.append(td);
    });
    tr.addEventListener("click", () => {
      [...table.querySelectorAll("tr")].forEach((el) => el.classList.remove("on"));
      tr.classList.add("on");
      showSite(row);
    });
    table.append(tr);
  });
  if (rows[0]) showSite(rows[0]);
  limits.replaceChildren();
  (payload.limitations || []).forEach((item) => {
    const li = document.createElement("li");
    li.textContent = item;
    limits.append(li);
  });
}

async function runAudit() {
  runButton.disabled = true;
  trace.replaceChildren();
  traceMode.textContent = "Live audit trace.";
  const body = {
    gene: gene.value.trim(),
    population: population.value,
    note: note.value,
    guides: guides.map((row, i) => ({
      id: row.idInput.value.trim() || `g${i + 1}`,
      spacer: row.spacerInput.value.trim(),
    })).filter((guide) => guide.spacer),
  };
  const response = await fetch("/api/audit", {
    method: "POST",
    headers: { "content-type": "application/json" },
    body: JSON.stringify(body),
  });
  const reader = response.body.getReader();
  const decoder = new TextDecoder();
  let buffer = "";
  while (true) {
    const chunk = await reader.read();
    if (chunk.done) break;
    buffer += decoder.decode(chunk.value, { stream: true });
    const parts = buffer.split("\n\n");
    buffer = parts.pop();
    parts.forEach((part) => {
      const line = part.replace(/^data: /, "");
      if (!line.trim()) return;
      const event = JSON.parse(line);
      addTrace(event);
      if (event.stage === "result") renderResult(event.payload);
    });
  }
  runButton.disabled = false;
  refreshSavedRuns().catch((error) => {
    traceMode.textContent = `Audit finished; saved-run list could not refresh: ${error}`;
  });
}

async function replaySavedRun() {
  if (!savedRun.value) return;
  loadRunButton.disabled = true;
  trace.replaceChildren();
  traceMode.textContent = "Loading saved trace from SQLite; no external API calls are made.";
  try {
    const response = await fetch(`/api/runs/${encodeURIComponent(savedRun.value)}`);
    if (!response.ok) throw new Error(`Could not load saved run (${response.status}).`);
    const data = await response.json();
    (data.events || []).forEach((event) => {
      addTrace(event);
      if (event.stage === "result") renderResult(event.payload);
    });
    if (!(data.events || []).some((event) => event.stage === "result")) {
      renderResult(data.record || {});
    }
    traceMode.textContent = `Offline replay · ${data.run.gene || "Unknown"} / ${data.run.population || "unknown"} · stored ${new Date(data.run.created_at).toLocaleString()}`;
  } catch (error) {
    traceMode.textContent = `Saved trace could not be replayed: ${error}`;
  } finally {
    loadRunButton.disabled = false;
  }
}

runButton.addEventListener("click", () => runAudit().catch((error) => {
  addTrace({ stage: "revise", title: "Browser error", detail: String(error) });
  runButton.disabled = false;
}));

loadRunButton.addEventListener("click", () => replaySavedRun());

document.querySelector("#refuse").addEventListener("click", () => {
  note.value = "Please diagnose this patient's genome and tell me if it is safe to inject.";
  runAudit();
});

document.querySelector("#bench").addEventListener("click", async () => {
  benchOut.textContent = "Running the CPU scan against chr11…";
  const response = await fetch("/api/bench", { method: "POST" });
  const data = await response.json();
  const lines = [
    data.gpu_note || data.error,
    data.cpu_seconds != null ? `CPU PAM-skip: ${data.cpu_seconds.toFixed(3)} s, ${Math.round(data.bases_per_sec).toLocaleString()} bases/s, ${data.hits} hits` : "",
    data.scan_seconds != null ? `CPU full scan: ${data.scan_seconds.toFixed(3)} s` : "",
    data.algorithms_equal != null ? `Hit lists equal: ${data.algorithms_equal}` : "",
    data.gpu && data.gpu.available ? `GPU: ${data.gpu.seconds.toFixed(3)} s, equal=${data.gpu.equal}` : "GPU: not timed on this machine.",
  ];
  benchOut.textContent = lines.filter(Boolean).join("\n");
});
