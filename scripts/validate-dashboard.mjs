import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { fileURLToPath } from 'node:url';
import path from 'node:path';
import vm from 'node:vm';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const read = file => readFileSync(path.join(root, file), 'utf8');
const rows = JSON.parse(read('data/literature/eukaryotes.snapshot.json'));
const provenance = JSON.parse(read('data/literature/provenance.json'));
assert.equal(rows.length, provenance.records);
const elements = new Map();
const element = id => {
  if (!elements.has(id)) {
    const classes = new Set();
    elements.set(id, { value: id === 'mapColor' ? 'SampleType' : id === 'mapSymbol' ? 'MolecularMethod' : '', style: {}, textContent: '', innerHTML: '', checked: false, dataset: {}, children: [], classList: { add: x => classes.add(x), remove: x => classes.delete(x), contains: x => classes.has(x) }, querySelectorAll() { return []; }, appendChild(node) { this.children.push(node); }, replaceChildren() { this.children = []; }, setAttribute() {}, addEventListener() {} });
  }
  return elements.get(id);
};
const plots = new Map();
const context = vm.createContext({
  console, Set, Map, URL, Blob, AbortController, setTimeout, clearTimeout,
  window: { HOLI_SNAPSHOT: rows, addEventListener() {} },
  document: { getElementById: element, querySelectorAll: () => [], createElement: () => element(Symbol()), head: { appendChild() {} } },
  Plotly: { setPlotConfig() {}, react(el, data, layout) { el.classList.add('js-plotly-plot'); plots.set(el, { data, layout }); }, purge(el) { plots.delete(el); el.classList.remove('js-plotly-plot'); } },
  localStorage: { getItem: () => null }, alert() {}, fetch: async () => ({ ok: false }),
});
vm.runInContext(read('dashboard/data-utils.js'), context);
vm.runInContext(read('dashboard/app.js'), context);
const run = code => vm.runInContext(code, context);
assert.equal(element('sSites').textContent, rows.length);
assert.equal(Number(element('sPubs').textContent), new Set(rows.map(row => context.HoliData.pubKey(row))).size);
assert.equal(plots.size, 5, 'All five existing plot families render');
const mapped = plots.get(element('map')).data.reduce((sum, trace) => sum + trace.lat.length, 0);
assert.equal(mapped, rows.filter(context.HoliData.validCoordinates).length);

const D = context.HoliData;
const parsed = D.parseDelimited('\uFEFFDOI,SiteName,year_published,Publication\r\n10.1/a,"Lake, A",2020,"Line 1\nLine 2"');
assert.equal(parsed.length, 1); assert.equal(parsed[0].SiteName, 'Lake, A'); assert.equal(parsed[0].Publication, 'Line 1\nLine 2');
assert.equal(D.parseDelimited('DOI\tSiteName\tyear_published\n10.1/a\tLake\t2020')[0].SiteName, 'Lake');
assert.throws(() => D.parseDelimited('Issue,Listed DOI\nNeeds review,10.1/a'), /Required columns/);
assert.throws(() => D.parseDelimited('DOI,SiteName,year_published\n10.1/a,"Bad,2020'), /Unclosed/);
assert.equal(D.validCoordinates({ Latitude: '', Longitude: '0' }), false);
assert.equal(D.validCoordinates({ Latitude: '91', Longitude: '0' }), false);
assert.equal(D.validCoordinates({ Latitude: '1N', Longitude: '0' }), false);
assert.equal(D.validCoordinates({ Latitude: '0', Longitude: '0' }), true);

run(`ROWS = [
 {DOI:'10.1234/ABC',Publication:'Paper A',year_published:'2020',SiteName:'Lake A',MolecularMethod:'SG',reference_target:'Ref',Latitude:'50',Longitude:'10'},
 {DOI:'https://doi.org/10.1234/abc',Publication:'Paper A',year_published:'2020',SiteName:'Lake B',MolecularMethod:'SG',reference_target:'Ref',Latitude:'95',Longitude:'10'},
 {DOI:'10.1234/DEF',Publication:'Paper B',year_published:'2021',SiteName:'Lake C',MolecularMethod:'TE',Latitude:'',Longitude:''}
]; initFromData(); render();`);
assert.equal(Number(element('sPubs').textContent), 2, 'Normalised DOI variants represent one paper');
assert.equal(Number(element('sSites').textContent), 3, 'Unmapped records remain in totals');
assert.equal(plots.get(element('map')).data.reduce((sum, trace) => sum + trace.lat.length, 0), 1);
assert.equal(run('countByYearMethod(ROWS).SG[2020]'), 1, 'Multi-site records deduplicate within year/method');
assert.equal(run("countDistinct(ROWS,'reference_target').Ref"), 1);
run(`searchTerm='no match'; render();`);
assert.equal(Number(element('sSites').textContent), 0);
assert.match(element('tableSummary').textContent, /No source records/);
run(`searchTerm=''; ROWS=[]; initFromData(); render();`);
assert.equal(run('state.yMin'), 0, 'Empty input does not create infinite year bounds');
run(`currentSource='Local file';`);
await run('refreshFromSheet(true)');
assert.match(element('srcLabel').innerHTML, /Local file/);
assert.doesNotMatch(element('srcLabel').innerHTML, /embedded snapshot/i, 'Failed refresh retains actual source identity');

const html = read('index.html');
assert.ok(html.startsWith('<!doctype html>'), 'No stray Markdown fences');
const ids = [...html.matchAll(/\bid="([^"]+)"/g)].map(m => m[1]);
assert.equal(new Set(ids).size, ids.length, 'Unique element IDs');
for (const match of html.matchAll(/(?:src|href)="([^"#]+)"/g)) {
  if (/^(https?:|mailto:)/.test(match[1])) continue;
  assert.ok(existsSync(path.join(root, match[1].split('#')[0])), `Missing asset: ${match[1]}`);
}
for (const name of ['SG_TE_map_method.png','SG_TE_map_SampleType_Target.png','barplot_no_publications.png','barplot_cumsum_no_publications_methods.png','barplot_no_publications_methods.png','barplot_N_databases.png','barplot_N_mappers.png','Holi_Workflow2025.png']) {
  const original = execFileSync('git', ['show', `main:figures/${name}`], { cwd: root, maxBuffer: 10e6 });
  const current = readFileSync(path.join(root, 'figures', name));
  assert.equal(createHash('sha256').update(original).digest('hex'), createHash('sha256').update(current).digest('hex'), `Original figure changed: ${name}`);
}
const originalHTML = execFileSync('git', ['show', 'main:index.html'], { cwd: root, encoding: 'utf8', maxBuffer: 10e6 });
const snapshotMatch = originalHTML.match(/const ROWS_EMBED = (\[.*?\]);/s);
if (snapshotMatch) assert.deepEqual(rows, JSON.parse(snapshotMatch[1]), 'Snapshot preserved from main');
console.log(`Passed: ${rows.length} source rows, 5 plots, DOI deduplication, coordinates, CSV parsing, empty results, failed-refresh provenance, assets, and all 8 original figures.`);
