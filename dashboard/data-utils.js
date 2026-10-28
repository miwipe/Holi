/* Shared data rules for the browser and dependency-free validation. */
(function (root) {
  'use strict';
  const clean = value => value == null ? '' : String(value).trim();
  const normalizeDOI = value => clean(value).replace(/^https?:\/\/(?:dx\.)?doi\.org\//i, '').replace(/^doi:\s*/i, '').toLowerCase();
  const pubKey = row => normalizeDOI(row.DOI) || clean(row.Publication).toLowerCase();
  const number = value => clean(value) !== '' && Number.isFinite(Number(value)) ? Number(value) : null;
  const validCoordinates = row => {
    const lat = number(row.Latitude), lon = number(row.Longitude);
    return lat !== null && lon !== null && Math.abs(lat) <= 90 && Math.abs(lon) <= 180;
  };
  function parseDelimited(text) {
    text = text.replace(/^\uFEFF/, '').replace(/\r\n/g, '\n').replace(/\r/g, '\n');
    const delim = text.split('\n')[0].includes('\t') ? '\t' : ',';
    const rows = []; let field = '', row = [], quoted = false;
    for (let i = 0; i < text.length; i++) {
      const ch = text[i];
      if (quoted) {
        if (ch === '"' && text[i + 1] === '"') { field += '"'; i++; }
        else if (ch === '"') quoted = false;
        else field += ch;
      } else if (ch === '"') quoted = true;
      else if (ch === delim) { row.push(field); field = ''; }
      else if (ch === '\n') { row.push(field); rows.push(row); row = []; field = ''; }
      else field += ch;
    }
    if (quoted) throw new Error('Unclosed quoted field in the table.');
    if (field || row.length) { row.push(field); rows.push(row); }
    const header = (rows.shift() || []).map(clean);
    if (!header.includes('DOI') || !header.includes('SiteName') || !header.includes('year_published')) throw new Error('Required columns: DOI, SiteName, year_published. Load the Eukaryotes table, not the DOI review tab.');
    return rows.filter(values => values.some(v => clean(v))).map(values => Object.fromEntries(header.map((key, i) => [key, clean(values[i])])));
  }
  root.HoliData = { clean, normalizeDOI, pubKey, number, validCoordinates, parseDelimited };
})(globalThis);
