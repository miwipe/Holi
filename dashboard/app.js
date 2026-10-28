
const ROWS_EMBED = window.HOLI_SNAPSHOT;

const PALETTE = ["#E69F00","#56B4E9","#009E73","#8A2BE2","#0072B2","#D55E00","#CC79A7",
  "#999999","#A52A2A","#00CED1","#FFD700","#228B22","#DA70D6","#006D6F","#6A5ACD","#C44E52","#556B2F"];
const SYMBOLS = ["circle","diamond","square","triangle-up","cross","star","pentagon","hexagon","x"];
const T = { ice:"#173538", frost:"#385555", dim:"#536c6c", line:"#e0e6dd", core:"#f5f5ef" };
const FONT = '"Inter", sans-serif';

let ROWS = ROWS_EMBED.slice();
let colorMap = {};                 // value -> color, stable across charts
function colorFor(v){ if(!(v in colorMap)){ colorMap[v]=PALETTE[Object.keys(colorMap).length % PALETTE.length]; } return colorMap[v]; }

const num = HoliData.number;
const clean = HoliData.clean;
const escapeHTML = value => clean(value).replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
let searchTerm = '', currentPage = 0;
const pageSize = 25;
let currentSource = 'Embedded snapshot';


// ---------- filter state ----------
const state = { method:new Set(), group:new Set(), sample:new Set(), yMin:0, yMax:0 };

function distinct(field){ return [...new Set(ROWS.map(r=>clean(r[field])).filter(Boolean))].sort(); }
function yearVals(){ return ROWS.map(r=>parseInt(r.year_published)).filter(y=>!isNaN(y)); }

function buildChips(){
  const mk=(host,vals,key,rustOn)=>{
    host.innerHTML='';
    vals.forEach(v=>{
      const c=document.createElement('button');
      c.type='button'; c.setAttribute('aria-pressed','true');
      c.className='chip on'+(rustOn?' rust':''); c.textContent=v;
      // chips start selected; clicking moves the value in/out of the excluded set
      c.onclick=()=>{ if(state[key].has(v)){ state[key].delete(v); c.classList.add('on'); }
                      else { state[key].add(v); c.classList.remove('on'); } c.setAttribute('aria-pressed',String(!state[key].has(v))); render(); };
      host.appendChild(c);
    });
  };
  mk(document.getElementById('fMethod'), distinct('MolecularMethod'),'method',true);
  mk(document.getElementById('fGroup'),  distinct('TargetGroup'),'group',false);
  mk(document.getElementById('fSample'), distinct('SampleType'),'sample',false);
}
document.querySelectorAll('[data-all]').forEach(b=>{
  b.onclick=()=>{ const k=b.dataset.all;
    const host={method:'fMethod',group:'fGroup',sample:'fSample'}[k];
    state[k].clear();                                   // exclude nothing = select all
    document.querySelectorAll('#'+host+' .chip').forEach(c=>{c.classList.add('on');c.setAttribute('aria-pressed','true');});
    render();
  };
});

function setupYears(){
  const ys=yearVals(); const lo=ys.length?Math.min(...ys):0, hi=ys.length?Math.max(...ys):0;
  const yMin=document.getElementById('yMin'), yMax=document.getElementById('yMax');
  [yMin,yMax].forEach(s=>{s.min=lo;s.max=hi;s.step=1;});
  yMin.value=lo; yMax.value=hi; state.yMin=lo; state.yMax=hi;
  const sync=()=>{ let a=+yMin.value,b=+yMax.value; if(a>b){[a,b]=[b,a];}
    state.yMin=a;state.yMax=b;
    document.getElementById('yMinLbl').textContent=a;
    document.getElementById('yMaxLbl').textContent=b; render(); };
  yMin.oninput=sync; yMax.oninput=sync;
  document.getElementById('yMinLbl').textContent=lo;
  document.getElementById('yMaxLbl').textContent=hi;
}
document.getElementById('resetAll').onclick=()=>{ searchTerm=''; document.getElementById('searchInput').value=''; initFromData(); render(); };
document.getElementById('searchInput').oninput=e=>{searchTerm=e.target.value.trim().toLowerCase();render();};

// filter state holds the DE-SELECTED (excluded) values per facet;
// empty set = nothing excluded = everything shown (including blank cells)
function passFilters(r){
  if(searchTerm && !['Publication','SiteName','RegionCountry','DOI','SampleType'].some(k=>clean(r[k]).toLowerCase().includes(searchTerm))) return false;
  const y=parseInt(r.year_published);
  if(!isNaN(y) && (y<state.yMin||y>state.yMax)) return false;
  if(state.method.has(clean(r.MolecularMethod))) return false;
  if(state.group.has(clean(r.TargetGroup))) return false;
  if(state.sample.has(clean(r.SampleType))) return false;
  return true;
}

// ---------- aggregation (dedupe pubs by DOI) ----------
const pubKey = HoliData.pubKey;
function countDistinct(rows, groupField){      // distinct (group, pub) -> count per group
  const seen=new Set(), out={};
  rows.forEach(r=>{ let g=clean(r[groupField]); if(!g) g='(unspecified)';
    const k=g+'||'+pubKey(r); if(seen.has(k))return; seen.add(k);
    out[g]=(out[g]||0)+1; });
  return out;
}
function countByYearMethod(rows){
  const seen=new Set(), out={};
  rows.forEach(r=>{ const y=parseInt(r.year_published); if(isNaN(y))return;
    const m=clean(r.MolecularMethod)||'?'; const k=y+'|'+m+'|'+pubKey(r);
    if(seen.has(k))return; seen.add(k);
    (out[m]=out[m]||{})[y]=(out[m][y]||0)+1; });
  return out;
}

// ---------- layout helpers ----------
function baseLayout(extra){
  return Object.assign({
    paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)',
    font:{color:T.frost,family:FONT,size:12},
    margin:{l:54,r:16,t:34,b:56},
    legend:{font:{color:T.frost,size:11},bgcolor:'rgba(0,0,0,0)'},
    xaxis:{gridcolor:T.line,zerolinecolor:T.line,linecolor:T.line,tickfont:{color:T.dim}},
    yaxis:{gridcolor:T.line,zerolinecolor:T.line,linecolor:T.line,tickfont:{color:T.dim}},
    title:{font:{color:T.ice,family:'Georgia,serif',size:15},x:0,xanchor:'left'}
  }, extra||{});
}
const CONF_BASE={responsive:true,displaylogo:false,displayModeBar:true,
  modeBarButtonsToRemove:['lasso2d','select2d','zoomIn2d','zoomOut2d','autoScale2d','toggleSpikelines'],
  toImageButtonOptions:{format:'png',scale:2}};
function confFor(name){
  return Object.assign({},CONF_BASE,
    {toImageButtonOptions:Object.assign({},CONF_BASE.toImageButtonOptions,{filename:name})});
}
// Only clear the container when it's NOT already a Plotly graph (e.g. it holds the
// "no data" placeholder). Wiping a live graph's DOM and then calling Plotly.react on it
// can throw — especially for the scattergeo map. When a graph exists, let react update it.
function prepPlot(el){
  if(!el) return;
  el.querySelectorAll(':scope > .empty').forEach(node=>node.remove());
  if(!el.classList.contains('js-plotly-plot')) el.innerHTML='';
}

// ---------- MAP ----------
function drawMap(rows){
  const el=document.getElementById('map');
  const colorBy=document.getElementById('mapColor').value;
  const symBy=document.getElementById('mapSymbol').value;
  const labels=document.getElementById('mapLabels').checked;
  const countAll=rows.length; rows=rows.filter(HoliData.validCoordinates);
  document.getElementById('mapCoverage').textContent=rows.length+' of '+countAll+' filtered site records have usable coordinates.';
  if(!rows.length){ Plotly.purge(el); el.innerHTML='<div class="empty">No mapped sites match the current filter.</div>'; return; }
  prepPlot(el);
  const symVals = symBy ? [...new Set(rows.map(r=>clean(r[symBy])))].sort() : [];
  const symOf = v => symVals.length? SYMBOLS[symVals.indexOf(v)%SYMBOLS.length] : 'circle';
  const groups={};
  rows.forEach(r=>{ const g=clean(r[colorBy])||'(unspecified)'; (groups[g]=groups[g]||[]).push(r); });
  const traces=Object.keys(groups).sort().map(g=>{
    const rs=groups[g];
    return {
      type:'scattergeo', mode: labels?'markers+text':'markers', name:g,
      lon:rs.map(r=>num(r.Longitude)), lat:rs.map(r=>num(r.Latitude)),
      text:rs.map(r=>clean(r.SiteName)),
      textposition:'top center', textfont:{size:8,color:T.frost},
      marker:{ size:8, color:colorFor(g), opacity:.9, line:{width:.6,color:'#ffffff'},
               symbol: symBy? rs.map(r=>symOf(clean(r[symBy]))) : 'circle' },
      hovertemplate:rs.map(r=>
        '<b>'+escapeHTML(r.SiteName)+'</b><br>'+escapeHTML(r.RegionCountry)+
        '<br>'+escapeHTML(r.Publication)+
        '<br>'+escapeHTML(r.SampleType)+' · '+escapeHTML(r.MolecularMethod)+
        ((clean(r.OldestYear_BP)&&clean(r.OldestYear_BP)!=='0')?('<br>oldest: '+escapeHTML(r.OldestYear_BP)+' yr BP'):'')+
        '<extra></extra>')
    };
  });
  const layout=baseLayout({
    margin:{l:0,r:0,t:6,b:0}, height:520,
    geo:{ projection:{type:'robinson'}, bgcolor:'rgba(0,0,0,0)',
      showland:true, landcolor:'#dce5d8', showcountries:true, countrycolor:'#b4c4b5',
      coastlinecolor:'#a8bcae', showcoastlines:true, showframe:false,
      showocean:true, oceancolor:'#f0f5f2', lataxis:{showgrid:true,gridcolor:'#dee8e0'},
      lonaxis:{showgrid:true,gridcolor:'#dee8e0'} },
    legend:{font:{color:T.frost,size:11},bgcolor:'rgba(255,255,255,.88)',y:.5}
  });
  Plotly.react(el,traces,layout,confFor('holi_map'));
}

// ---------- YEAR (stacked by method) ----------
function drawYear(rows){
  const el=document.getElementById('cYear');
  const by=countByYearMethod(rows);
  const ys=yearVals(); const lo=state.yMin, hi=state.yMax;
  const years=[]; for(let y=lo;y<=hi;y++) years.push(y);
  const methods=Object.keys(by).sort();
  if(!methods.length){ Plotly.purge(el); el.innerHTML='<div class="empty">No publications match.</div>'; return; }
  prepPlot(el);
  const traces=methods.map(m=>({
    type:'bar', name:m, x:years, y:years.map(y=>(by[m][y]||0)),
    marker:{color:colorFor(m)} }));
  const layout=baseLayout({ barmode:'stack', title:{text:''},
    xaxis:{gridcolor:T.line,linecolor:T.line,tickfont:{color:T.dim},dtick:1,title:{text:'Year',font:{color:T.dim,size:11}}},
    yaxis:{gridcolor:T.line,linecolor:T.line,tickfont:{color:T.dim},title:{text:'Publications',font:{color:T.dim,size:11}}},
    legend:{orientation:'h',y:1.12,x:0} });
  Plotly.react(el,traces,layout,confFor('holi_publications_by_year_method'));
}

// ---------- CUMULATIVE ----------
function drawCum(rows){
  const el=document.getElementById('cCum');
  const by=countByYearMethod(rows);
  const lo=state.yMin, hi=state.yMax; const years=[]; for(let y=lo;y<=hi;y++) years.push(y);
  const methods=Object.keys(by).sort();
  if(!methods.length){ Plotly.purge(el); el.innerHTML='<div class="empty">No publications match.</div>'; return; }
  prepPlot(el);
  const traces=methods.map(m=>{ let c=0; const cum=years.map(y=>{c+=(by[m][y]||0);return c;});
    return {type:'scatter',mode:'lines',name:m,x:years,y:cum,stackgroup:'one',
      line:{width:.5,color:colorFor(m)},fillcolor:colorFor(m)+'cc'}; });
  const layout=baseLayout({ title:{text:''},
    xaxis:{gridcolor:T.line,linecolor:T.line,tickfont:{color:T.dim},dtick:1,title:{text:'Year',font:{color:T.dim,size:11}}},
    yaxis:{gridcolor:T.line,linecolor:T.line,tickfont:{color:T.dim},title:{text:'Cumulative',font:{color:T.dim,size:11}}},
    legend:{orientation:'h',y:1.12,x:0} });
  Plotly.react(el,traces,layout,confFor('holi_cumulative_by_method'));
}

// ---------- horizontal count bars ----------
function drawCountBar(elId,rows,field,title,color){
  const el=document.getElementById(elId);
  const counts=countDistinct(rows,field);
  let entries=Object.entries(counts).sort((a,b)=>a[1]-b[1]);
  if(!entries.length){ Plotly.purge(el); el.innerHTML='<div class="empty">No data in current filter.</div>'; return; }
  prepPlot(el);
  el.style.height=Math.max(390,entries.length*22+85)+'px';
  const trace={ type:'bar', orientation:'h',
    y:entries.map(e=>e[0]), x:entries.map(e=>e[1]),
    marker:{color:color}, text:entries.map(e=>e[1]), textposition:'outside',
    textfont:{color:T.frost,size:11},
    hovertemplate:'%{y}<br>%{x} publications<extra></extra>' };
  const layout=baseLayout({ title:{text:''},
    margin:{l:170,r:30,t:34,b:40},
    xaxis:{gridcolor:T.line,linecolor:T.line,tickfont:{color:T.dim},title:{text:'Publications',font:{color:T.dim,size:11}}},
    yaxis:{linecolor:T.line,tickfont:{color:T.frost,size:11},automargin:true} });
  Plotly.react(el,[trace],layout,confFor('holi_'+title.toLowerCase().replace(/[^a-z0-9]+/g,'_').replace(/^_|_$/g,'')));
}

// ---------- stat strip ----------
function updateStats(rows){
  const pubs=new Set(rows.map(pubKey).filter(Boolean));
  const ys=rows.map(r=>parseInt(r.year_published)).filter(y=>!isNaN(y));
  const ages=rows.map(r=>num(r.OldestYear_BP)).filter(v=>v&&v>0);
  document.getElementById('sPubs').textContent=pubs.size;
  document.getElementById('sSites').textContent=rows.length;
  document.getElementById('sYears').textContent = ys.length? (Math.min(...ys)+'–'+Math.max(...ys)) : '–';
  document.getElementById('sOldest').textContent = ages.length? Math.max(...ages).toLocaleString() : '–';
}

// ---------- core-log decoration ----------
function paintCore(){
  const el=document.getElementById('coreCol'); if(!el) return; el.innerHTML='';
  const bands=[['#c3b087',26],['#a78b60',16],['#7c9b95',30],['#4b716e',12],['#6e8a72',20],
               ['#2a3947',24],['#5a6b3f',14],['#33424f',28],['#22303a',18]];
  bands.forEach(([c,h])=>{ const b=document.createElement('i'); b.style.height=h+'px';
    b.style.background=c; b.style.boxShadow='inset 0 1px 0 rgba(255,255,255,.05)'; el.appendChild(b); });
}

// ---------- master render ----------
function render(){
  const rows=ROWS.filter(passFilters);
  updateStats(rows);
  currentPage=0; drawRecords(rows);
  if(typeof Plotly==='undefined'){ showPlotlyBanner(); return; }
  const safe=(fn,name)=>{ try{ fn(); }catch(e){ console.error('render('+name+') failed:',e); } };
  safe(()=>drawMap(rows),'map');
  safe(()=>drawYear(rows),'year');
  safe(()=>drawCum(rows),'cum');
  safe(()=>drawCountBar('cDb',rows,'reference_target','Publications by reference database','#a78bdb'),'db');
  safe(()=>drawCountBar('cMap',rows,'mapper','Publications by mapper','#5fb3a3'),'mapper');
}

// ---------- init / data load ----------
function initFromData(){
  colorMap={};
  ['MolecularMethod','SampleType','TargetGroup','Lab'].forEach(field=>distinct(field).forEach(colorFor));
  state.method=new Set(); state.group=new Set(); state.sample=new Set();
  buildChips(); setupYears();
  document.getElementById('srcCount').textContent=ROWS.length;
}

const parseDelimited=HoliData.parseDelimited;
document.getElementById('fileIn').onchange=e=>{
  const f=e.target.files[0]; if(!f)return;
  const rd=new FileReader();
  rd.onload=()=>{ try{
      const parsed=parseDelimited(rd.result);
      const ok=parsed;
      if(!ok.length){ alert('No source records found. Check the headers.'); return; }
      ROWS=ok;
      currentSource='Local file'; document.getElementById('sourceBadge').textContent=currentSource;
      setSrc('Source: '+escapeHTML(f.name)+' · <b id="srcCount">'+ROWS.length+'</b> site records');
      initFromData(); render();
    }catch(err){ alert('Could not parse file: '+err.message); }
  };
  rd.readAsText(f);
};

// ---------- live Google Sheet ----------
const SHEET_ID = "13cmBUi4cigUaTKtQeFLFvS0gXT8AeWxWKzHv2UcOBCI";
let   SHEET_GID = "0";   // gid of the DATA tab — the number after #gid= in the sheet's URL
// "Publish to web → CSV" link (CORS-friendly). This is what loads live by default.
const PUBLISHED_CSV_URL = "https://docs.google.com/spreadsheets/d/e/2PACX-1vQXmAeJPrjrUQtRSqeipWDUA5zanFtrb_tOwkyzTHm2TkWCv493GH4wzInB-E-Ny5XxiUxYCmesAh8g/pub?output=csv";
const SHEET_URL_KEY = "holi_sheet_csv_url";

// restore a previously pasted "Publish to web" CSV link, if any
try{ const saved=localStorage.getItem(SHEET_URL_KEY);
     if(saved){ const f=document.getElementById('sheetUrl'); if(f) f.value=saved; } }catch(e){}
(function(){ const f=document.getElementById('sheetUrl');
  if(f) f.addEventListener('change',()=>{ try{ localStorage.setItem(SHEET_URL_KEY, f.value.trim()); }catch(e){} }); })();

function sheetUrls(){
  const urls=[];
  const f=document.getElementById('sheetUrl');
  const pasted=f && f.value.trim();
  if(pasted) urls.push(pasted);                                  // optional override
  if(PUBLISHED_CSV_URL) urls.push(PUBLISHED_CSV_URL);            // the published default (CORS-friendly)
  urls.push(`https://docs.google.com/spreadsheets/d/${SHEET_ID}/gviz/tq?tqx=out:csv&gid=${SHEET_GID}`);
  urls.push(`https://docs.google.com/spreadsheets/d/${SHEET_ID}/export?format=csv&gid=${SHEET_GID}`);
  return urls;
}
function setSrc(html){ const el=document.getElementById('srcLabel'); if(el) el.innerHTML=html; }
function sourceFailed(note){ setSrc('Source: '+escapeHTML(currentSource)+' · <b id="srcCount">'+ROWS.length+'</b> site records · '+escapeHTML(note)); }
async function fetchSheetRows(){
  for(const url of sheetUrls()){
    try{
      const ctrl=new AbortController(); const t=setTimeout(()=>ctrl.abort(),9000);
      let res; try { res=await fetch(url,{signal:ctrl.signal,credentials:'omit'}); } finally { clearTimeout(t); }
      if(!res.ok) continue;
      const text=await res.text();
      if(/^\s*</.test(text)) continue;              // an HTML login/error page, not CSV
      const parsed=parseDelimited(text);
      const ok=parsed;
      if(ok.length) return ok;
    }catch(e){ /* CORS / offline / aborted — try next endpoint */ }
  }
  return null;
}
async function refreshFromSheet(silent){
  setSrc('Checking Google Sheet… · keeping <b id="srcCount">'+ROWS.length+'</b> '+escapeHTML(currentSource)+' records');
  const rows=await fetchSheetRows();
  if(rows){
    ROWS=rows; currentSource='Live Google Sheet'; document.getElementById('sourceBadge').textContent=currentSource; initFromData(); render();
    setSrc('Source: live Google Sheet · <b id="srcCount">'+ROWS.length+'</b> site records · loaded '+new Date().toLocaleString());
  } else {
    sourceFailed('live refresh unavailable; previous data retained');
    if(!silent) alert(
      'Could not load the live Google Sheet.\n\n'+
      'A browser usually cannot read a sheet just from "Anyone with the link" — Google blocks it (CORS).\n\n'+
      'To enable live loading:\n'+
      '  1. Open the sheet → File → Share → Publish to web\n'+
      '  2. Choose the data tab and the "Comma-separated values (.csv)" format\n'+
      '  3. Copy the link it gives you (it contains "/pub?...output=csv")\n'+
      '  4. Paste it into the box next to the Refresh button\n\n'+
      'Keeping the previously loaded dataset.');
  }
}
document.getElementById('sheetRefresh').onclick=()=>refreshFromSheet(false);

// ---------- download the (filtered) literature table ----------
function downloadFilteredCSV(){
  const rows=ROWS.filter(passFilters);
  if(!rows.length){ alert('No rows match the current filters — nothing to download.'); return; }
  // union of columns across rows, preserving first-seen order
  const seen=new Set(), cols=[];
  rows.forEach(r=>Object.keys(r).forEach(k=>{ if(!seen.has(k)){ seen.add(k); cols.push(k);} }));
  const esc=v=>{ v=(v==null?'':String(v)); return /[",\n]/.test(v)?'"'+v.replace(/"/g,'""')+'"':v; };
  const lines=[cols.map(esc).join(',')];
  rows.forEach(r=>lines.push(cols.map(c=>esc(r[c])).join(',')));
  const blob=new Blob(['\ufeff'+lines.join('\r\n')],{type:'text/csv;charset=utf-8'});
  const url=URL.createObjectURL(blob);
  const a=document.createElement('a');
  a.href=url; a.download='holi_filtered_'+rows.length+'_sites.csv';
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(()=>URL.revokeObjectURL(url),1500);
}
document.getElementById('tableDownload').onclick=downloadFilteredCSV;
document.getElementById('recordsDownload').onclick=downloadFilteredCSV;

// ---------- table and explicit exports ----------
function drawRecords(rows){
  const pages=Math.ceil(rows.length/pageSize);
  currentPage=Math.max(0,Math.min(currentPage,Math.max(0,pages-1)));
  const host=document.getElementById('recordRows');host.replaceChildren();
  rows.slice(currentPage*pageSize,(currentPage+1)*pageSize).forEach(row=>{
    const tr=document.createElement('tr');
    ['Publication','SiteName','RegionCountry','MolecularMethod','SampleType','DOI'].forEach(field=>{
      const td=document.createElement('td');const value=clean(row[field]);
      if(field==='DOI' && /^10\.\d{4,9}\/\S+$/.test(HoliData.normalizeDOI(value))){
        const a=document.createElement('a');a.href='https://doi.org/'+encodeURI(HoliData.normalizeDOI(value));a.textContent=value;a.target='_blank';a.rel='noopener';td.appendChild(a);
      }else td.textContent=value || '—';
      tr.appendChild(td);
    });host.appendChild(tr);
  });
  document.getElementById('tableSummary').textContent=rows.length ? `${currentPage*pageSize+1}–${Math.min((currentPage+1)*pageSize,rows.length)} of ${rows.length} site records` : 'No source records match the filters.';
  document.getElementById('prevPage').disabled=currentPage===0;
  document.getElementById('nextPage').disabled=currentPage+1>=pages;
}
document.getElementById('prevPage').onclick=()=>{currentPage--;drawRecords(ROWS.filter(passFilters));};
document.getElementById('nextPage').onclick=()=>{currentPage++;drawRecords(ROWS.filter(passFilters));};
document.querySelectorAll('[data-export]').forEach(button=>{
  button.onclick=async()=>{
    const el=document.getElementById(button.dataset.export);
    if(typeof Plotly==='undefined'||!el.classList.contains('js-plotly-plot')) {alert('No chart available to export for the current filter.');return;}
    button.disabled=true;
    try{await Plotly.downloadImage(el,{format:'png',filename:'holi_'+button.dataset.export,width:1200,height:button.dataset.export==='map'?700:650,scale:2});}
    catch(error){alert('PNG export failed: '+error.message);}
    finally{button.disabled=false;}
  };
});
function showPlotlyBanner(){
  if(document.getElementById('plotlyBanner'))return;
  const p=document.createElement('p');p.id='plotlyBanner';p.className='review-notice';p.textContent='The charting library could not be loaded. Source records and filters remain available. Check that dashboard/vendor is included with this page.';
  document.getElementById('explore').prepend(p);
}
function loadPlotly(onReady){
  if(typeof Plotly==='undefined'){showPlotlyBanner();return;}
  Plotly.setPlotConfig({topojsonURL:'dashboard/vendor/'});
  onReady();
}

// boot — UI scaffolding first (works without Plotly), then charts once it loads
paintCore();
initFromData();
updateStats(ROWS.filter(passFilters));
['mapColor','mapSymbol','mapLabels'].forEach(id=>
  document.getElementById(id).addEventListener('change',()=>{
    if(typeof Plotly!=='undefined') drawMap(ROWS.filter(passFilters)); }));
window.addEventListener('resize',()=>{ /* plotly responsive handles it */ });
loadPlotly(()=>{
  render();                                  // embedded snapshot renders immediately
  const auto=document.getElementById('sheetAuto');
  if(auto && auto.checked) refreshFromSheet(true);   // then quietly try the live sheet
});
