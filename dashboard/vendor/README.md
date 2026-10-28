# Dashboard vendor assets

- `plotly-2.35.2.min.js`: Plotly.js v2.35.2, extracted from the previous repository `index.html`. Copyright 2012–2024 Plotly, Inc.; distributed under the MIT licence. The upstream distribution and notices are retained. See `plotly-LICENSE.txt`.
- `world_110m.json`: map topology fetched on 27 October 2026 from https://cdn.plot.ly/world_110m.json, the asset used by Plotly's geographic traces. The basemap is derived from Natural Earth public-domain data. See https://www.naturalearthdata.com/about/terms-of-use/.

These are pinned local assets, so displaying charts and the default world map does not require a live CDN. No vendor code was rewritten as part of the redesign.
