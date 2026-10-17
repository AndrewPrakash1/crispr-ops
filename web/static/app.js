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

const guides = [];

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
}

runButton.addEventListener("click", () => runAudit().catch((error) => {
  addTrace({ stage: "revise", title: "Browser error", detail: String(error) });
  runButton.disabled = false;
}));

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
