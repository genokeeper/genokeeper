"use strict";
// The dashboard's script. It looks the chosen view up in the page's data and builds the
// feed from it with createElement and textContent; it never filters by level or topic,
// because `genokeeper.dashboard` embedded one query result per threshold pair and topic.
// The search box only narrows a view to the items whose embedded `search` string holds it.
// The Changes tab looks each change's kind up the same way, per threshold pair and base, and
// narrows the table by its chips. Every text it shows is a message of `data.messages`, the
// catalog `resources/dashboard/messages.json` in the page's language.

// Every table here and in the data has no prototype, so no lookup by a data string or a hash
// value (a topic, a cause, `#constructor=1`) finds an inherited name.
const bare = (value) => (value && typeof value === "object" && !Array.isArray(value) ? Object.assign(Object.create(null), value) : value);
const data = JSON.parse(document.getElementById("data").textContent, (key, value) => bare(value));
const M = data.messages;
const CARDINAL = new Intl.PluralRules(data.lang);
const NUMBER = new Intl.NumberFormat(data.lang);
const number = (n) => NUMBER.format(n);
const ORDINAL = new Intl.PluralRules(data.lang, { type: "ordinal" });

// The message `key` with its `{name}` placeholders filled from `params`, in one pass. A
// message with forms takes the one `rules` select for `params.n`, else its `other` (the
// catalog's Russian keeps its "many" form there).
function t(key, params = {}, rules = CARDINAL) {
  const message = M[key];
  const text = typeof message === "string" ? message : message[rules.select(params.n)] ?? message.other;
  return text.replace(/\{(\w+)\}/g, (_, name) => String(params[name]));
}

const SCALES = {
  reliability: ["low", "medium", "high"],
  validity: ["conflicting", "limited", "moderate", "strong"],
};
// Each level's name and its short form under a bar.
const LEVELS = M.levels;
// A caller's banner is the demo's frame (`genokeeper.demo`), the one caller that sets one.
const DEMO = document.querySelector(".frame-note") !== null;
const TOPICS = M.topics;
// A repeat count's class (`genokeeper.repeats.loci.RepeatClass`).
const REPEAT_CLASSES = M.repeat_classes;

// The page's one URL: a PubMed link is this prefix and a PMID of digits only, which
// `genokeeper.history.read_publications` checked.
const PUBMED = "https://pubmed.ncbi.nlm.nih.gov/";

const [reliability, validity, topic] = data.default.split("|");
// `chosen`: a hash or a tier chip set the Changes tab's thresholds or tier.
const state = { reliability, validity, topic, search: "", tab: "feed", chosen: false };

function el(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

function labelled(label, ...rest) {
  const span = el("span");
  span.append(el("b", "", `${label}: `), ...rest);
  return span;
}

function contig(f) {
  return f.marker.grch38 === null ? null : f.marker.grch38.split(":")[0];
}

// A score stored since slice 3f shows its decile (one genome's score is less precise than
// its percentile); an older one its percentile.
function percentile(effect) {
  const score = (effect || {}).score;
  if (score && score.decile !== undefined && score.decile !== null) return t("decile", { n: score.decile }, ORDINAL);
  return score && score.percentile !== null ? t("percentile", { n: score.percentile }, ORDINAL) : t("no_percentile");
}

// Slice 3f's trait statements (`genokeeper.sources.traits`) and their results.
const TRAITS = M.traits;
const TRAIT_RESULTS = M.trait_results;

function percent(share) {
  return `${Math.round(100 * share)}%`;
}

// ClinVar's strings read as words, `_` a space and `|` a comma; the data stays as stored,
// and another source's text (`pgsc_calc`, a transcript id) is its own.
function readable(text, source = "clinvar") {
  return typeof text === "string" && source === "clinvar" ? text.replace(/_/g, " ").replace(/\|/g, ", ") : text;
}

function joined(...parts) {
  return parts.filter((part) => part !== null && part !== undefined && part !== "").join(", ");
}

// ClinVar's condition names that name no condition, as read (`readable`) and lower-cased.
const PLACEHOLDERS = new Set(["not provided", "not specified", "see cases", "none provided", "not applicable"]);

// A ClinVar claim's subtitle: a drug-response row's drugs (its conditions `<drug>_response…`)
// and what grades it (A3s: never its stars); else its classification and first three named conditions.
function clinvarHeadline(f) {
  const { classification, condition } = f.claim;
  const parts = String(condition || "").split("|");
  const drugs = [...new Set(parts.map((part) => /^(.+?)_response(?:_-_|$)/i.exec(part)?.[1]).filter(Boolean)
    .map((drug) => readable(drug).toLowerCase()))];
  if (classification === "drug_response" && drugs.length) {
    return [t("drug_response", { drugs: drugs.join(", ") }), f.validity.reason].filter(Boolean).join("; ");
  }
  const named = parts.map((part) => readable(part)).filter((part) => part && !PLACEHOLDERS.has(part.toLowerCase()));
  const more = named.length > 3 ? t("more", { n: named.length - 3 }) : "";
  return joined(readable(classification), named.slice(0, 3).join(", ")) + more;
}

// A slice 3f claim's result (HLA, MITOMAP, HHV-6, ROH, blood groups, ancestry, traits), or
// undefined for any other claim type.
function typingHeadline(claim) {
  switch (claim.type) {
    case "hla_allele":
      return joined(t("hla_allele", claim), claim.condition);
    case "hla_association":
      return joined(t("hla_association", claim), claim.condition);
    case "hla_coeliac": {
      const t1k = claim.t1k || {};
      const types = [["DQ2.5", t1k.dq25], ["DQ2.2", t1k.dq22], ["DQ8", t1k.dq8]].filter(([, n]) => n).map(([name]) => name);
      const found = types.length ? t("coeliac_found", { types: types.join(", ") }) : t1k.half ? t("coeliac_half") : t("coeliac_none");
      return t("coeliac", { found });
    }
    case "hla_table": {
      return t("hla_table", { n: (claim.genes || []).length });
    }
    case "mitomap":
      return joined(`MITOMAP ${claim.allele}`, claim.condition, typeof claim.level === "number" ? t("mitomap_level", { level: percent(claim.level) }) : null);
    case "hhv6_ici":
      return joined(t("hhv6", { species: claim.species || "HHV-6" }), typeof claim.ratio === "number" ? t("hhv6_ratio", claim) : null);
    case "roh_summary": {
      return t("roh", { n: claim.n_long || 0, mb: ((claim.total_long || 0) / 1e6).toFixed(2) });
    }
    case "blood_group":
      return t("blood_group", { system: claim.system, value: claim.phenotype || claim.genotype || t("unknown_f") });
    case "ancestry_similarity":
      return t("ancestry_similarity", { population: claim.population || t("unknown_f") }) + (claim.low_confidence ? t("low_confidence") : "");
    case "ancestry_admixture": {
      const shares = Object.entries(claim.proportions || {}).sort(([, a], [, b]) => b - a);
      return t("ancestry_admixture", { shares: shares.map(([group, share]) => `${group} ${percent(share)}`).join(", ") || t("unknown_pl") });
    }
    case "ancestry_neanderthal":
      return typeof claim.share === "number" ? t("neanderthal", { share: (100 * claim.share).toFixed(2) }) : t("no_neanderthal");
    case "trait": {
      const result = claim.result;
      const shown = TRAIT_RESULTS[result] ?? (typeof result === "number" ? t("effect_copies", { n: result }) : result ?? t("unknown_m"));
      return `${TRAITS[claim.trait] || claim.trait}: ${shown}`;
    }
    default:
      return undefined;
  }
}

// The row's result: each claim type's own, else its gene call's, else its condition.
function headline(f) {
  const claim = f.claim || {};
  const call = f.gene_call;
  const typed = typingHeadline(claim);
  if (typed !== undefined) return typed;
  if (claim.type === "clinvar_germline") return clinvarHeadline(f);
  // ClinPGx's level is the grade (`sources.pharmgkb` writes "ClinPGx level <level>").
  const clinpgx = /^ClinPGx level (\S+)$/.exec(f.validity.source_level || "");
  if (claim.type === "pharmgkb_clinical_annotation" && clinpgx && (claim.drugs || []).length) {
    return t("drug_response_clinpgx", { drugs: claim.drugs.map((d) => d.drug).join(", "), level: clinpgx[1] });
  }
  if (claim.type === "mito_haplogroup") return t("mito_haplogroup", { haplogroup: claim.haplogroup || t("unknown_f") });
  if (claim.type === "y_haplogroup") return t("y_haplogroup", { haplogroup: claim.haplogroup || t("unknown_f") });
  if (claim.type === "pgs") return `${claim.condition}: ${percentile(f.effect)}${claim.weak ? `, ${t("weak_power")}` : ""}`;
  if (claim.type === "apoe_diplotype") return t("apoe", { diplotype: claim.diplotype || t("unknown_m") });
  if (claim.type === "predicted_variant") return joined(claim.classification, claim.consequence, claim.hgvsp || claim.hgvsc);
  if (claim.type === "repeat_expansion" && call) {
    return joined(call.diplotype || t("no_size"), REPEAT_CLASSES[claim.class] || claim.class, call.phenotype);
  }
  if (call) return `${call.diplotype || t("no_diplotype")}, ${call.phenotype || t("no_phenotype")}`;
  return claim.condition || claim.type;
}

function markerLine(f) {
  const { marker, genotype } = f;
  const parts = [];
  if (marker.grch38 === null) {
    parts.push(t("gene_level"));
  } else {
    parts.push(marker.rsid || `GRCh38 ${marker.grch38}`);
    if (contig(f) === "chrM") {
      // mutserve's dosage 1 is presence in the mitochondria, at the level its warning gives.
      parts.push(t(genotype.dosage === "0" ? "chrM_absent" : "chrM_present", genotype));
      parts.push(...f.reliability.warnings.filter((w) => w.startsWith("heteroplasmy")));
    } else {
      parts.push(t("copies", genotype));
    }
  }
  parts.push(readable(f.validity.source_level, f.sources[0].source));
  return parts.filter((part) => part).join(" · ");
}

function bar(label, levels, level, kind) {
  const track = el("span", `bar ${kind}`);
  const fill = el("i");
  fill.style.width = `${(100 * (levels.indexOf(level) + 1)) / levels.length}%`;
  track.append(fill);
  return [el("span", "", label), track, el("span", "", LEVELS[level][1])];
}

function genotypeText(f) {
  const call = f.gene_call;
  if (call) {
    const callers = Object.entries(call.callers || {}).map(([name, c]) => `${name} ${c}`);
    return [
      t("reliability_is", { level: LEVELS[call.level][0] }),
      t("callers", { callers: callers.join("; ") || t("none") }),
      ...call.warnings,
    ].join(", ");
  }
  const r = f.reliability;
  if (r.level === null) return t("no_genotype");
  let pWrong = t("not_measured");
  if (r.p_wrong !== null) {
    const ci = r.p_wrong_ci ? ` [${r.p_wrong_ci.map((p) => p.toExponential(1)).join(", ")}]` : "";
    pWrong = t("p_wrong_measured", { p: r.p_wrong.toExponential(1), ci });
  }
  return [
    f.genotype.gt,
    t("reliability_is", { level: LEVELS[r.level][0] }),
    t("p_wrong", { p: pWrong }),
    r.gq === null ? null : `GQ ${r.gq}`,
    r.row_type,
    t("flags", { flags: r.flags.join(", ") || t("none") }),
    ...r.warnings,
  ].filter((part) => part !== null).join(", ");
}

function effectText(effect) {
  const parts = [];
  const { relative, absolute, score } = effect;
  if (relative) {
    const ci = relative.ci ? ` [${relative.ci.join(", ")}]` : "";
    parts.push(`${relative.measure} ${relative.value}${ci}`);
  }
  if (absolute) {
    parts.push(t("risk", absolute));
  }
  if (score) {
    const z = (score.z === null ? "" : `, z ${score.z}`) + (score.z_norm2 === undefined || score.z_norm2 === null ? "" : `, Z_norm2 ${score.z_norm2}`);
    parts.push(t("score", { ...score, percentile: percentile(effect), z }));
  }
  return parts.join("; ");
}

// A REViewer image is a file on this machine. Its anchor comes from the `file://` literal
// here, which no stored text can move because the CSP pins this script by hash and the
// path only follows the prefix, percent-encoded; a path that is not absolute gets no link.
function imageLink(claim) {
  if (typeof claim.image !== "string" || !claim.image.startsWith("/")) return null;
  const a = el("a", "", t("reads_image"));
  a.href = `file://${encodeURI(claim.image)}`;
  a.target = "_blank";
  a.rel = "noopener noreferrer";
  return labelled(t("reads"), a);
}

function opened(f) {
  const box = el("div", "open");
  box.append(el("b", "", t("what_it_means")));
  const e = f.explanation;
  if (e) {
    // "claude-opus-5-5" reads "Opus 5.5", "claude-sonnet-4-20250514" "Sonnet 4"; another id stays as stored.
    const id = /^claude-([a-z]+)-(\d+)(?:-(\d{1,2}))?(?:-\d{8})?$/.exec(e.model);
    const model = id ? `${id[1][0].toUpperCase()}${id[1].slice(1)} ${id[2]}${id[3] ? `.${id[3]}` : ""}` : e.model;
    const tip = t("explanation_tooltip", { model, demo: DEMO ? t("for_this_demo") : "" });
    // An engaging entry's gist leads; every generated text carries the tooltip.
    const parts = [...(e.gist ? [el("span", "gist", e.gist)] : []), el("span", "", e.text),
      el("span", "note", t("written_by", { model: e.model, cites: e.cites.join(", ") }))];
    for (const part of parts) {
      part.setAttribute("title", tip);
      part.setAttribute("aria-description", tip);
    }
    box.append(...parts);
  } else {
    box.append(el("span", "", t("no_explanation")));
  }
  box.append(labelled(t("genotype"), genotypeText(f)));
  const v = f.validity;
  const own = v.source_level ? ` (${readable(v.source_level, f.sources[0].source)})` : "";
  const reason = v.reason ? `; ${v.reason}` : "";
  box.append(labelled(t("evidence"), `${LEVELS[v.level][0]}${own}${reason}`));
  if (f.effect) box.append(labelled(t("effect"), effectText(f.effect)));
  const image = imageLink(f.claim || {});
  if (image) box.append(image);
  box.append(labelled(t("actionable"), t(f.actionable ? "yes" : "no")));
  if (f.dates) box.append(labelled(t("dates"), datesText(f.dates)));
  const sources = el("ul");
  for (const s of f.sources) sources.append(el("li", "", t("source_version", s)));
  box.append(el("b", "", t("sources")), sources);
  return box;
}

function row(f) {
  const details = el("details", "row");
  const summary = el("summary");
  const main = el("div");
  // A result with no gene (a score, a haplogroup, ancestry) reads as its own headline.
  const raw = f.gene ? `${f.gene}: ${headline(f)}` : headline(f);
  const e = f.explanation;
  // The model's plain headline leads; the source's own words stay just under it.
  // The model sometimes names the gene in its title already.
  const named = !f.gene || e?.title.toUpperCase().includes(f.gene.toUpperCase());
  if (e) main.append(el("b", "", named ? e.title : `${f.gene}: ${e.title}`), el("div", "note", raw));
  else main.append(el("b", "", raw));
  main.append(el("div", "sub", markerLine(f)));
  if (["chrX", "chrY"].includes(contig(f))) main.append(el("div", "warn", t("haploid_note")));
  const pair = el("div", "pair");
  const level = f.reliability.level || (f.gene_call && f.gene_call.level);
  if (level) pair.append(...bar(t("bar_genotype"), SCALES.reliability, level, "rel"));
  pair.append(...bar(t("bar_evidence"), SCALES.validity, f.validity.level, "ev"));
  summary.append(main, pair);
  details.append(summary, opened(f));
  return details;
}

function reasonText(reason) {
  const unassessed = reason.startsWith("unassessed:");
  const [name, threshold] = reason.replace(/^unassessed:/, "").split(":");
  const level = LEVELS[threshold] ? LEVELS[threshold][0] : threshold;
  const key = bare({
    reliability_below: "reliability_below",
    validity_below: "validity_below",
    topic: "other_topics",
    research_association: "research_hidden",
  })[name];
  if (key === undefined) return reason;
  const text = t(key, { level });
  return unassessed ? t("unknown_genotype", { text }) : text;
}

// The pair fewest steps from the current one whose view of `topic` shows rows (the looser on
// a tie): its thresholds, its count and the labels of the thresholds that differ; or null.
function nearest(topic) {
  const steps = (r, v) => Math.abs(SCALES.reliability.indexOf(r) - SCALES.reliability.indexOf(state.reliability)) +
    Math.abs(SCALES.validity.indexOf(v) - SCALES.validity.indexOf(state.validity));
  const pairs = SCALES.reliability.flatMap((r) => SCALES.validity.map((v) => ({ r, v, n: count(`${r}|${v}|${topic}`) })));
  const best = pairs.filter((p) => p.n).sort((a, b) => steps(a.r, a.v) - steps(b.r, b.v))[0];
  if (!best) return null;
  const where = joined(best.r === state.reliability ? null : M.matrix_rows[best.r], best.v === state.validity ? null : M.matrix_columns[best.v]);
  return { ...best, n: number(best.n), where };
}

// `ids` by their topic, in the catalog's order of topics, each topic once.
function byTopic(ids) {
  const groups = new Map(Object.keys(TOPICS).map((name) => [name, []]));
  for (const id of ids) {
    const name = data.items[id].topic;
    if (!groups.has(name)) groups.set(name, []);
    groups.get(name).push(id);
  }
  return [...groups].filter(([, list]) => list.length);
}

function chips() {
  const present = new Set(Object.values(data.items).map((f) => f.topic));
  const known = Object.keys(TOPICS).filter((name) => present.has(name));
  const box = document.getElementById("chips");
  box.replaceChildren();
  for (const name of ["", ...known, ...[...present].filter((name) => !(name in TOPICS))]) {
    const label = name ? TOPICS[name] || name : t("all_topics");
    // A topic with no row at this cell says where its nearest rows are, and moves there.
    const near = count(`${state.reliability}|${state.validity}|${name}`) ? null : nearest(name);
    const chip = el("button", name === state.topic ? "chip on" : "chip", near ? t("chip_empty", { topic: label, ...near }) : label);
    chip.type = "button";
    chip.addEventListener("click", () => {
      Object.assign(state, { topic: name }, near && { reliability: near.r, validity: near.v });
      render();
      writeHash();
    });
    box.append(chip);
  }
}

// The run's flagged non-findings (`genokeeper.dashboard._excluded`): a screening result
// that rules something out, each with the test that would confirm it.
function renderExcluded() {
  const box = document.getElementById("excluded");
  const list = el("ul");
  for (const x of data.excluded) list.append(el("li", "", t("excluded_item", { ...x, gene: x.gene || t("no_gene") })));
  if (!data.excluded.length) list.append(el("li", "", t("no_exclusions")));
  box.replaceChildren(el("b", "", t("excluded")), list);
}

// The run's other stage-3 watch items and non-findings (`genokeeper.dashboard._screened`),
// which the feed never holds: a watch item is neither normal nor a finding, so each kind is
// its own list.
function renderScreened() {
  const box = document.getElementById("screened");
  box.replaceChildren();
  for (const [kind, heading] of [["watch", "screened_watch"], ["normal", "screened_normal"]]) {
    const list = el("ul");
    for (const x of data.screened[kind]) {
      const item = el("li", "", x.wording === null ? x.label : t("screened_item", x));
      // A watch item's explanation, written as one (`genokeeper.explain.WATCH_NOTE`).
      const e = x.explanation;
      if (e) item.append(el("div", "", e.text), el("span", "note", t("written_by", { model: e.model, cites: e.cites.join(", ") })));
      list.append(item);
    }
    if (!data.screened[kind].length) list.append(el("li", "", t("none_in_run")));
    box.append(el("b", "", t(heading)), list);
  }
}

// A view's findings that hold the typed text. Each item's `search` is
// `genokeeper.query.searchable`'s, and so is the typed text. The rows are searched once per
// typed text, and each view's count kept until the text changes.
let matched = { search: null, ids: null, counts: new Map() };
function searched(view) {
  if (matched.search !== state.search) {
    const q = state.search.toLowerCase().replaceAll("_", " ").replaceAll("|", ", ");
    const ids = new Set(Object.keys(data.items).filter((id) => data.items[id].search.includes(q)));
    matched = { search: state.search, ids, counts: new Map() };
  }
  return view.findings.filter((id) => matched.ids.has(id));
}

// The number of findings the view `key` (`<reliability>|<validity>|<topic>`) shows.
function count(key) {
  if (matched.search !== state.search || !matched.counts.has(key)) {
    const n = searched(data.views[key]).length;
    matched.counts.set(key, n);
  }
  return matched.counts.get(key);
}

// What a reader would act on or ask about: ClinVar and predicted findings, a non-typical drug
// or carrier result, a score in the top or bottom decile; never a trait or ancestry (the
// haplogroups' topic).
const HANDFUL = 5;
function standsOut(f) {
  const claim = f.claim || {};
  if (f.normal || f.kind !== "finding" || ["traits", "ancestry"].includes(f.topic)) return false;
  if (claim.type === "pgs") return [1, 10].includes(f.effect?.score?.decile);
  return ["clinvar_germline", "predicted_variant"].includes(claim.type) || ["drugs", "carrier"].includes(f.topic);
}

function renderFeed() {
  const view = data.views[`${state.reliability}|${state.validity}|${state.topic}`];
  chips();

  const shown = searched(view);
  const total = Object.keys(data.items).length;
  document.getElementById("count").textContent = t("count", { shown: number(shown.length), n: total, total: number(total) });
  // Typical results (`normal`) fold into one line per topic, open while a search runs; the
  // others are named by topic above the feed, each opening its row.
  const rows = new Map(shown.map((id) => [id, row(data.items[id])]));
  const standing = shown.filter((id) => !data.items[id].normal);
  const folds = byTopic(shown.filter((id) => data.items[id].normal)).map(([name, ids]) => {
    const fold = el("details", "fold");
    fold.open = Boolean(state.search);
    fold.append(el("summary", "", t("typical", { n: ids.length, topic: TOPICS[name] || name })), ...ids.map((id) => rows.get(id)));
    return fold;
  });
  document.getElementById("feed").replaceChildren(...standing.map((id) => rows.get(id)), ...folds);
  const out = document.getElementById("stands-out");
  const outstanding = standing.filter((id) => standsOut(data.items[id]));
  out.replaceChildren(...(outstanding.length ? [el("b", "", t("stands_out"))] : []));
  for (const [name, ids] of byTopic(outstanding)) {
    const line = el("div", "", `${TOPICS[name] || name}: `);
    const labels = new Map();
    for (const id of ids) {
      const label = data.items[id].gene || headline(data.items[id]);
      if (!labels.has(label)) labels.set(label, id);
    }
    for (const [i, [label, id]] of [...labels].slice(0, HANDFUL).entries()) {
      const link = el("button", "link", label);
      link.type = "button";
      link.addEventListener("click", () => {
        rows.get(id).open = true;
        rows.get(id).scrollIntoView?.();
      });
      line.append(...(i ? [", "] : []), link);
    }
    if (labels.size > HANDFUL) {
      const more = el("button", "link", t("more", { n: labels.size - HANDFUL }));
      more.type = "button";
      more.addEventListener("click", () => {
        state.topic = name;
        render();
        writeHash();
      });
      line.append(more);
    }
    out.append(line);
  }
  renderExcluded();
  renderScreened();

  // One clause per reason, the research layer's with or without a genotype as one.
  const hidden = new Map();
  for (const [reason, n] of Object.entries(view.hidden)) {
    const text = reason.endsWith("research_association") ? reasonText("research_association") : reasonText(reason);
    hidden.set(text, (hidden.get(text) || 0) + n);
  }
  const unmatched = view.findings.length - shown.length;
  if (unmatched) hidden.set(t("unmatched"), unmatched);
  const note = document.getElementById("hidden");
  note.replaceChildren();
  if (hidden.size) {
    const count = number([...hidden.values()].reduce((sum, n) => sum + n, 0));
    const list = el("ul");
    for (const [text, n] of hidden) list.append(el("li", "", `${number(n)} ${text}`));
    note.append(t("hidden", { count, how: t(unmatched ? "pick_cell_or_search" : "pick_cell") }), list);
  }

  const unassessed = document.getElementById("unassessed");
  unassessed.replaceChildren();
  if (view.unassessed.total) {
    const list = el("ul");
    for (const [gene, entry] of Object.entries(view.unassessed.by_gene)) {
      const reasons = Object.entries(entry.reasons).map(([reason, n]) => `${reason} (${n})`);
      list.append(el("li", "", `${gene || t("no_gene")}: ${entry.sites}; ${reasons.join("; ")}`));
    }
    unassessed.append(
      el("b", "", t("unassessed", { n: view.unassessed.total })),
      list,
    );
  }

  // What the run never looked at (`genokeeper.carrier`'s lists), so no absence reads as
  // "not a carrier".
  const uncovered = document.getElementById("not-covered");
  uncovered.replaceChildren();
  if (state.topic === "carrier") {
    const all = el("ul");
    for (const text of data.not_covered.all) all.append(el("li", "", text));
    const genes = el("ul");
    for (const [gene, text] of Object.entries(data.not_covered.genes)) genes.append(el("li", "", `${gene}: ${text}`));
    uncovered.append(
      el("b", "", t("not_checked")),
      all,
      el("b", "", t("uncovered_genes")),
      genes,
    );
  }
}

// --- the Changes tab ------------------------------------------------------------------------

const changes = data.changes;
// Each change kind's cell, chip and cell class, in the order the table sorts them.
const KINDS = bare(Object.fromEntries(
  Object.entries({ entered: "up", left: "down", level: "mid", new: "mid", gone: "mid", details: "mid", publications: "mid" })
    .map(([kind, className]) => [kind, [...M.kinds[kind], className]]),
));
// Run inputs (in a cause), refreshed sources (in the banner) and the update check's
// sources, by their names.
const SOURCES = bare({
  clinvar38: "ClinVar",
  clinvar: "ClinVar",
  clinvar_citations: t("clinvar_citations"),
  gene_validity: "ClinGen",
  clingen: "ClinGen",
  clingen_dosage: "ClinGen dosage",
  gencc: "GenCC",
  cgd: "CGD",
  faf: "gnomAD FAF",
  "pgx.cpic": "CPIC",
  cpic: "CPIC",
  "pgx.pharmgkb": "ClinPGx",
  pharmgkb: "ClinPGx",
  "pgx.positions": "PharmCAT",
  gwas: "GWAS Catalog",
  vep: "VEP",
  medlineplus: "MedlinePlus",
  litvar: "LitVar",
  pharmcat: "PharmCAT",
  ipd_imgt_hla: "IPD-IMGT/HLA",
  panelapp: "PanelApp",
  ensembl_vep: "Ensembl VEP",
  haplogrep3: "Haplogrep3",
  pgs_catalog: "PGS Catalog",
  gnomad: "gnomAD",
});
const STATEMENT_KINDS = M.statement_kinds;
const FIELDS = M.fields;
const tab = {
  base: changes.default_base,
  tier: "first",
  kind: "",
  topic: "",
  cause: "",
  sort: null,
  descending: false,
  open: new Set(),
};

function level(name) {
  return name === null || name === undefined ? "—" : LEVELS[name] ? LEVELS[name][0] : name;
}

// A date on the oldest snapshot, for a statement already in it, is "on or before" it.
function day(dates, name, before) {
  const date = dates[name];
  if (date === null) return "—";
  return dates.on_or_before && date === changes.oldest ? `${before} ${date}` : date;
}

function datesText(dates) {
  return [
    t("found_on", { date: day(dates, "first_found", t("on_or_before")) }),
    t("changed_on", { date: day(dates, "last_changed", t("on_or_before")) }),
    t("assessment_changed_on", { date: day(dates, "assessment_changed", t("on_or_before")) }),
  ].join(", ");
}

function causeGroup(cause) {
  const own = bare({ analyser: "cause_analyser", genotype: "cause_genotype", reference: "cause_reference" })[cause];
  if (own) return t(own);
  const input = cause.replace(/^source:/, "").split("@")[0];
  return SOURCES[input] || input;
}

function causeText(cause) {
  const release = cause.startsWith("source:") ? cause.split("@").slice(1).join("@") : "";
  return release ? `${causeGroup(cause)} ${release}` : causeGroup(cause);
}

function validityText(side, source) {
  if (side === null) return "—";
  const own = side.source_level ? ` (${readable(side.source_level, source)})` : "";
  return `${level(side.validity)}${own}`;
}

function reliabilityText(row) {
  const before = row.before && row.before.level;
  const after = row.after && row.after.level;
  if (row.before && row.after && before !== after) return `${level(before)} → ${level(after)}`;
  return level(after || before || null);
}

function marker(row) {
  return row.marker.rsid || row.marker.grch38 || t("gene_level");
}

function dated(row, name) {
  return row.dates ? day(row.dates, name, "≤") : "—";
}

// The table's columns: a title and a cell's text, which is also what a click sorts by.
const COLUMNS = [
  [M.columns.change, (row, kind) => KINDS[kind][0]],
  [M.columns.gene, (row) => `${row.gene || t("no_gene")} · ${marker(row)}`],
  [M.columns.topic, (row) => TOPICS[row.topic] || row.topic || "—"],
  [M.columns.validity, (row) => `${validityText(row.before, row.source)} → ${validityText(row.after, row.source)}`],
  [M.columns.reliability, reliabilityText],
  [M.columns.cause, (row) => row.causes.map(causeText).join("; ") || "—"],
  [M.columns.found, (row) => dated(row, "first_found")],
  [M.columns.changed, (row) => dated(row, "last_changed")],
  [M.columns.assessed, (row) => dated(row, "assessment_changed")],
];

function pubmed(pmid) {
  const a = el("a", "", `PMID ${pmid}`);
  a.href = PUBMED + pmid;
  a.target = "_blank";
  a.rel = "noopener noreferrer";
  return a;
}

function links(label, pmids) {
  const span = labelled(label);
  pmids.forEach((pmid, i) => span.append(...(i ? [", "] : []), pubmed(pmid)));
  return span;
}

function sideValue(side, name, source) {
  const value = side === null ? null : side[name];
  if (value === null || value === undefined) return "—";
  if (name === "validity" || name === "level") return level(value);
  if (name === "actionable") return t(value ? "yes" : "no");
  if (name === "kind") return STATEMENT_KINDS[value] || value;
  if (name === "source_level") return readable(value, source);
  return String(value);
}

function expanded(row) {
  const box = el("div", "open");
  if (row.fields.length) {
    const fields = el("ul");
    for (const name of row.fields) {
      const text = `${sideValue(row.before, name, row.source)} → ${sideValue(row.after, name, row.source)}`;
      fields.append(el("li", "", `${FIELDS[name] || name}: ${text}`));
    }
    box.append(el("b", "", t("what_changed")), fields);
  } else if (!row.before || !row.after) {
    box.append(el("span", "", t(row.before ? "statement_gone" : "statement_new")));
  }
  const timeline = el("ul");
  for (const h of row.history) {
    const parts = h.kind === null
      ? [t("none")]
      : [STATEMENT_KINDS[h.kind] || h.kind, t("history_evidence", { level: level(h.validity) }), t("history_genotype", { level: level(h.level) })];
    timeline.append(el("li", "", `${h.date} (${h.run}): ${parts.join(", ")}`));
  }
  box.append(el("b", "", t("history")), timeline);
  const before = row.before && row.before.source_version;
  const after = row.after && row.after.source_version;
  const version = before && after && before !== after ? `${before} → ${after}` : after || before || "—";
  box.append(labelled(t("source"), t("source_version", { source: row.source || "—", id: row.id, version })));
  if (row.causes.length) box.append(labelled(M.columns.cause, row.causes.map(causeText).join("; ")));
  if (row.publications.new.length) box.append(links(t("new_articles"), row.publications.new));
  if (row.publications.attached_by_lab.length) {
    box.append(links(t("attached_by_lab"), row.publications.attached_by_lab));
  }
  if (row.dates) box.append(labelled(t("dates"), datesText(row.dates)));
  return box;
}

function chipGroup(label, options, key) {
  const group = el("div", "group");
  group.append(el("span", "", label));
  for (const [value, text] of options) {
    const chip = el("button", tab[key] === value ? "chip on" : "chip", text);
    chip.type = "button";
    chip.addEventListener("click", () => {
      tab[key] = value;
      if (key === "tier") state.chosen = true;
      renderChanges();
    });
    group.append(chip);
  }
  return group;
}

function present(values, selected, order) {
  const found = new Set(values);
  if (selected) found.add(selected);
  return [...found].sort((a, b) => order(a) - order(b) || a.localeCompare(b, data.lang));
}

// A release that is a date (`20261108`) in ISO form; any other stays as stored.
function iso(release) {
  return /^\d{8}$/.test(release) ? `${release.slice(0, 4)}-${release.slice(4, 6)}-${release.slice(6)}` : release;
}

// The head: the two runs by their dates (their ids in tooltips) and the source versions
// that differ between them.
function renderHead(compared) {
  const head = document.getElementById("changes-head");
  const update = el("span", "", t("update_head", { date: changes.date || "—" }));
  update.setAttribute("title", changes.run);
  head.replaceChildren(update);
  if (!compared) return;
  const select = el("select");
  select.setAttribute("aria-label", t("compared_label"));
  select.setAttribute("title", tab.base);
  for (const [run, base] of Object.entries(changes.bases)) {
    const option = el("option", "", `${base.date} · ${t(base.refresh ? "base_refresh" : "base_run")}`);
    option.setAttribute("title", run);
    option.value = run;
    option.selected = run === tab.base;
    select.append(option);
  }
  select.addEventListener("change", () => {
    tab.base = select.value;
    tab.open.clear();
    renderChanges();
  });
  head.append(t("compared_with"), select);
  const differ = compared.sources.map((s) => `${SOURCES[s.name] || s.name} ${iso(s.before)} → ${iso(s.after)}`);
  if (differ.length) head.append(el("div", "", t("sources_differ", { sources: differ.join(", ") })));
}

// The update check's age is computed here, when the page is viewed: a check that stopped
// after the page was written still shows as stale. A demo lab has no check, and no block.
function renderCheck() {
  const box = document.getElementById("check");
  const info = changes.check;
  if (!info) {
    box.replaceChildren();
    return;
  }
  const ran = info.ran_at ? Date.parse(info.ran_at) : NaN;
  const days = Number.isNaN(ran) ? null : Math.max(0, Math.floor((Date.now() - ran) / 86400000));
  const notes = [];
  if (days === null || days > info.stale_days) {
    const warning = el("div", "banner");
    warning.append(
      el("b", "", days === null ? t("check_never") : t("check_stale", { n: days })),
      t("check_unknown"),
    );
    notes.push(warning);
  } else {
    const named = (items) => items.map((s) => `${SOURCES[s.name] || s.name} ${s.upstream}`).join(", ");
    const line = el("div", "note");
    line.textContent = t("check_ran", { n: days }) +
      (info.newer.length ? t("check_newer", { sources: named(info.newer) }) : t("check_none_newer"));
    notes.push(line);
  }
  if (info.failed.length && days !== null) {
    const unchecked = el("div", "note");
    unchecked.textContent = t("check_failed", { sources: info.failed.map((name) => SOURCES[name] || name).join(", ") });
    notes.push(unchecked);
  }
  box.replaceChildren(...notes);
}

function renderChanges() {
  const compared = changes.bases[tab.base];
  renderHead(compared);
  renderCheck();
  const banner = document.getElementById("banner");
  banner.replaceChildren(
    ...changes.failed.map((source) => {
      const box = el("div", "banner");
      const since = source.before ? t("since_old", { version: source.before }) : t("since_previous");
      box.append(
        el("b", "", t("not_updated", { source: SOURCES[source.name] || source.name })),
        t("not_updated_detail", { since }),
      );
      return box;
    }),
  );
  const note = document.getElementById("no-base");
  note.textContent = compared ? "" : t("no_base");
  const table = document.getElementById("changes");
  const chipsBox = document.getElementById("change-chips");
  const research = document.getElementById("research");
  table.hidden = chipsBox.hidden = !compared;
  research.textContent = "";
  if (!compared) return;

  const kinds = compared.views[`${state.reliability}|${state.validity}`];
  const kindOrder = (kind) => Object.keys(KINDS).indexOf(kind);
  const rows = Object.values(compared.rows);
  const first = rows.filter((row) => changes.first_tier.includes(kinds[row.id]));
  // Thresholds no one chose that keep no change, where the loosest keep some: all of them.
  const loosest = compared.views[`${SCALES.reliability[0]}|${SCALES.validity[0]}`];
  if (!state.chosen && !first.length && rows.some((row) => changes.first_tier.includes(loosest[row.id]))) {
    tab.tier = "all";
    note.textContent = t("below_thresholds", { n: rows.length });
  }
  const tiered = tab.tier === "first" ? first : rows;
  const shown = tiered.filter(
    (row) =>
      (!tab.kind || kinds[row.id] === tab.kind) &&
      (!tab.topic || row.topic === tab.topic) &&
      (!tab.cause || row.causes.some((cause) => causeGroup(cause) === tab.cause)),
  );
  const topics = Object.keys(TOPICS);
  chipsBox.replaceChildren(
    chipGroup(t("chip_level"), [["first", t("tier_first", { n: first.length })], ["all", t("tier_all", { n: rows.length })]], "tier"),
    chipGroup(M.columns.change, [["", t("chip_all")], ...present(tiered.map((row) => kinds[row.id]), tab.kind, kindOrder).map((k) => [k, KINDS[k][1]])], "kind"),
    chipGroup(M.columns.topic, [["", t("chip_all")], ...present(tiered.map((row) => row.topic || ""), tab.topic, (name) => topics.indexOf(name)).filter((name) => name).map((name) => [name, TOPICS[name] || name])], "topic"),
    chipGroup(M.columns.cause, [["", t("chip_all")], ...present(tiered.flatMap((row) => row.causes.map(causeGroup)), tab.cause, () => 0).map((c) => [c, c])], "cause"),
  );

  const cells = (row) => COLUMNS.map(([, text]) => text(row, kinds[row.id]));
  if (tab.sort === null) {
    shown.sort((a, b) => kindOrder(kinds[a.id]) - kindOrder(kinds[b.id]));
  } else {
    const sign = tab.descending ? -1 : 1;
    shown.sort((a, b) => sign * cells(a)[tab.sort].localeCompare(cells(b)[tab.sort], data.lang));
  }
  const header = el("tr");
  COLUMNS.forEach(([title], i) => {
    const arrow = tab.sort === i ? (tab.descending ? " ▼" : " ▲") : "";
    const th = el("th", "", title + arrow);
    th.addEventListener("click", () => {
      tab.descending = tab.sort === i && !tab.descending;
      tab.sort = i;
      renderChanges();
    });
    header.append(th);
  });
  const head = el("thead");
  head.append(header);
  const body = el("tbody");
  for (const row of shown) {
    const line = el("tr", "line");
    cells(row).forEach((text, i) => {
      const className = i === 0 ? KINDS[kinds[row.id]][2] : i >= 6 ? "date" : "";
      line.append(el("td", className, text));
    });
    line.addEventListener("click", () => {
      if (tab.open.has(row.id)) tab.open.delete(row.id);
      else tab.open.add(row.id);
      renderChanges();
    });
    body.append(line);
    if (tab.open.has(row.id)) {
      const more = el("tr");
      const cell = el("td");
      cell.colSpan = COLUMNS.length;
      cell.append(expanded(row));
      more.append(cell);
      body.append(more);
    }
  }
  if (!shown.length) {
    const none = el("tr");
    const cell = el("td", "note", t("no_changes"));
    cell.colSpan = COLUMNS.length;
    none.append(cell);
    body.append(none);
  }
  table.replaceChildren(head, body);

  if (compared.research.length) {
    const counts = compared.research.map((r) => {
      const why = r.causes.map(causeText).join("; ") || t("no_cause");
      return `${KINDS[r.kind][1]} — ${r.count} (${why})`;
    });
    research.textContent = t("research_counts", { counts: counts.join(", ") });
  }
}

function render() {
  for (const [cell, r, v] of CELLS) {
    const n = count(`${r}|${v}|${state.topic}`);
    const marked = cell.className ? [el("small", "", t("default_cell"))] : [];
    cell.replaceChildren(number(n), ...marked);
    cell.setAttribute("aria-pressed", String(r === state.reliability && v === state.validity));
    const name = t("cell_label", { n, count: number(n), row: M.matrix_rows[r], column: M.matrix_columns[v] });
    cell.setAttribute("aria-label", [name, ...marked.map((m) => m.textContent)].join(", "));
  }
  const onFeed = state.tab === "feed";
  for (const [name, on] of [["feed", onFeed], ["changes", !onFeed]]) {
    document.getElementById(`tab-${name}`).className = on ? "tab on" : "tab";
    document.getElementById(`${name}-view`).hidden = !on;
  }
  document.getElementById("search").hidden = document.getElementById("chips").hidden = !onFeed;
  if (onFeed) renderFeed();
  else renderChanges();
}

// The view in the URL's hash, `#tab=feed|changes&rel=0..2&ev=0..3&topic=<topic>&q=<text>`, so
// a link opens one exact view (from file:// too); an unknown key or value is ignored.
const HASH_SCALES = bare({ rel: "reliability", ev: "validity" });
function readHash() {
  for (const pair of location.hash.replace(/^#/, "").split("&")) {
    const at = pair.indexOf("=");
    if (at < 0) continue;
    const key = pair.slice(0, at);
    let value;
    try {
      value = decodeURIComponent(pair.slice(at + 1));
    } catch {
      continue; // a cut %-escape
    }
    if (key === "tab" && (value === "feed" || value === "changes")) state.tab = value;
    else if (key in HASH_SCALES && /^[0-9]$/.test(value)) {
      const levels = SCALES[HASH_SCALES[key]];
      if (Number(value) < levels.length) Object.assign(state, { [HASH_SCALES[key]]: levels[Number(value)], chosen: true });
    } else if (key === "topic" && `${state.reliability}|${state.validity}|${value}` in data.views) state.topic = value;
    else if (key === "q") state.search = value;
  }
}

// The view written back when a control moves, replacing the history entry; the default pair
// is left out, so a copied link opens the Changes tab as an unchosen one does.
function writeHash() {
  const parts = [`tab=${state.tab}`];
  if (state.reliability !== reliability || state.validity !== validity) {
    parts.push(`rel=${SCALES.reliability.indexOf(state.reliability)}`, `ev=${SCALES.validity.indexOf(state.validity)}`);
  }
  if (state.topic) parts.push(`topic=${encodeURIComponent(state.topic)}`);
  if (state.search) parts.push(`q=${encodeURIComponent(state.search)}`);
  history.replaceState(null, "", `#${parts.join("&")}`);
}

for (const name of ["feed", "changes"]) {
  document.getElementById(`tab-${name}`).addEventListener("click", () => {
    state.tab = name;
    render();
    writeHash();
  });
}

const search = document.getElementById("search");
state.search = search.value; // a value the browser kept across a reload
readHash();
search.value = state.search;
search.addEventListener("input", () => {
  state.search = search.value;
  render();
  writeHash();
});
// The matrix that sets both thresholds, as the demo's landing: genotype reliability down the
// side, "high only" on top, strength of evidence across; `render` counts each cell's view.
const matrix = document.getElementById("matrix");
matrix.setAttribute("aria-label", t("matrix_axis"));
// Each column head in full, and short for a phone's width (`dashboard.css`).
matrix.append(el("span", "", t("matrix_axis")), ...SCALES.validity.map((v) => {
  const head = el("span", "ev");
  head.append(el("span", "long", M.matrix_columns[v]), el("span", "short", M.matrix_columns_short[v]));
  return head;
}));
const CELLS = [...SCALES.reliability].reverse().flatMap((r) => {
  matrix.append(el("span", "rel", M.matrix_rows[r]));
  return SCALES.validity.map((v) => {
    const cell = el("button", r === reliability && v === validity ? "default" : "");
    cell.type = "button";
    cell.addEventListener("click", () => {
      Object.assign(state, { reliability: r, validity: v });
      render();
      writeHash();
    });
    matrix.append(cell);
    return [cell, r, v];
  });
});
// An edited hash is the view, read from the defaults as on load; `replaceState` fires no
// `hashchange`, so writing the view back never loops.
window.addEventListener("hashchange", () => {
  Object.assign(state, { reliability, validity, topic, search: "", tab: "feed", chosen: false });
  tab.tier = "first";
  readHash();
  search.value = state.search;
  render();
});
// The footer's first line (the run, ClinVar's release, genokeeper's version); every version
// and digest in a closed disclosure.
const { clinvar, genokeeper } = data.versions;
const first = [t("footer", { run: data.run_id }), clinvar && `ClinVar ${iso(clinvar)}`, genokeeper && `genokeeper ${genokeeper}`];
const versions = el("details");
versions.append(el("summary", "", t("versions_digests")), el("div", "", Object.entries(data.versions).map(([tool, version]) => `${tool} ${version}`).join(", ")));
document.getElementById("footer").replaceChildren(el("div", "", first.filter(Boolean).join(" · ")), versions);
render();
